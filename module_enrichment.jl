proj_info = (id = "vgirault_vzvapms",
             data_ver = "20180301",
             fit_ver = "20180301",
             model_obj = "protgroup",
             oesc_ver = "20201129")

using Pkg
Pkg.activate(@__DIR__)
#@everywhere (myid()==1) || (ENV["MKL_NUM_THREADS"]=1) # disable MKL threading on workers
using Distances, DataFrames, CategoricalArrays, CSV, RData, CodecZlib, JLD2

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

using OptEnrichedSetCover

include(joinpath(misc_scripts_path, "frame_utils.jl"));
include(joinpath(misc_scripts_path, "delimdata_utils.jl"));
include(joinpath(misc_scripts_path, "omics_collections.jl"));
include(joinpath(misc_scripts_path, "gmt_reader.jl"));

use_cached_jlser = true
if !use_cached_jlser

@info "Loading $scratch_path $(proj_info.fit_ver) MS GLM analysis data..."
msglm_rdata = load(joinpath(scratch_path, "$(proj_info.id)_msglm_fit_$(proj_info.fit_ver).RData"))
# option 1: filter all interactions directly
all_iactions_df = msglm_rdata["fit_contrasts"]["iactions"][:, :]
all_iactions_df = all_iactions_df[(all_iactions_df.var .== "iaction_labu_replCI") .&
                                  occursin.("_vs_allMinus", all_iactions_df.contrast) .&
                                  (all_iactions_df.condition_role .== "signal"), :]
#iactions_df = all_iactions_df[(all_iactions_df.prob_nonpos .<= 1E-2) .&
#                              (all_iactions_df.median_log2 .>= 0.5), :]
# option 2: use manually curated final list of interactions
iactions_report_df = CSV.read(joinpath(data_path, "$(proj_info.id)_edge-table_20191114.txt"), DataFrame)
iactions_df = iactions_report_df[:, :]
bait_col = :display_bait
categorical!(iactions_df, bait_col)

#include(joinpath(misc_scripts_path, "protgroup_utils.jl"));

const ProtgroupType = eltype(iactions_df.protgroup_id)

proteins_df = msglm_rdata["ms_data"]["proteins"][:, :];
protgroups_df = msglm_rdata["ms_data"]["protgroups"][:, :];
protgroups_df.protgroup_id = convert(Vector{ProtgroupType}, protgroups_df.protgroup_id)
protgroup2ac_df = msglm_rdata["ms_data"]["protein2protgroups"][:, :]
protgroup2ac_df.protgroup_id = convert(Vector{ProtgroupType}, protgroup2ac_df.protgroup_id)
protgroup2ac_noiso_df = protgroup2ac_df[:, :]
protgroup2ac_noiso_df.protein_ac = replace.(protgroup2ac_noiso_df.protein_ac, Ref(r"-\d+$" => ""))
protgroup2ac_df = unique!(vcat(protgroup2ac_df, protgroup2ac_noiso_df))
protgroup2ac_df = semijoin(protgroup2ac_df, proteins_df, on=:protein_ac) # only proteins from fasta (exclude contaminants and reverse)

#protgroup2ac_df = ProtgroupUtils.expand_dlm_column(protgroups_df);
prot_ac2pg_ids = Dict{String, Set{ProtgroupType}}();
for ac2ids_df in groupby(protgroup2ac_df[protgroup2ac_df.is_majority, :], :protein_ac)
    prot_ac2pg_ids[ac2ids_df.protein_ac[1]] = Set(ac2ids_df.protgroup_id)
end

# human mappings from http://download.baderlab.org/EM_Genesets/April_01_2019/Human/UniProt/
# FIXME using all evidence codes
genesets_df, genesets_coll = GMT.read(String,
        joinpath(party3rd_data_path, "Human_GO_AllPathways_with_GO_iea_October_01_2020_UniProt.gmt"),
        id_col = :term_id, src_col = :term_src);

pcomplexes_df, pcomplex_iactors_df, pcomplex_iactor2ac_df =
    OmicsCollections.ppicollection(joinpath(party3rd_data_path, "complexes_20191217.RData"), seqdb=:uniprot);
pcomplexes_df[!, :coll_id] .= "protein_complexes";

# make complexes collections, keep complexes with at least 2 participants
pcomplex_coll = FrameUtils.frame2collection(innerjoin(pcomplex_iactors_df, pcomplex_iactor2ac_df,
    on=[:file, :entry_index, :interaction_id, :interactor_id]),
            set_col=:complex_id, obj_col=:protein_ac, min_size=2)
protac_sets = merge(genesets_coll, pcomplex_coll)

terms_df = vcat(rename(genesets_df[!, [:term_src, :term_id, :name, :descr]],
                       :term_src => :coll_id, :name=>:term_name, :descr=>:term_descr),
                #rename(goterm_info_df[[:id, :name, :def]], :onto => :coll_id, :id=>:term_id, :name=>:term_name, :def=>:term_descr),
                rename(pcomplexes_df[!, [:coll_id, :complex_id, :interaction_label, :interaction_name]],
                       :complex_id=>:term_id, :interaction_label=>:term_name, :interaction_name=>:term_descr));
protac2term_df = FrameUtils.collection2frame(protac_sets, terms_df,
                                             setid_col=:term_id, objid_col=:protein_ac)

pg2term_df = select!(innerjoin(protgroup2ac_df, protac2term_df, on = :protein_ac),
                     Not([:protein_ac])) |> unique!
protac_colls = FrameUtils.frame2collections(protac2term_df, obj_col=:protein_ac,
                                            set_col=:term_id, coll_col=:coll_id)
pg_colls = FrameUtils.frame2collections(pg2term_df, obj_col=:protgroup_id,
                                        set_col=:term_id, coll_col=:coll_id)

using OBOParse

# GO ontology from http://geneontology.org/page/download-ontology
go_onto = OBOParse.load(joinpath(party3rd_data_path, "go_20201118.obo"), "GO");
goterm_info_df = OmicsCollections.terminfo_frame(go_onto)

proteincomplex_goid = "GO:0032991" # umbrella term for all protein complexes
# GO IDs for all protein complexes
proteincomplexes_goterms = descendants(go_onto, go_onto[proteincomplex_goid])

@save(joinpath(scratch_path, "$(proj_info.id)_omics_collections_$(proj_info.oesc_ver).jld2"),
      protgroups_df, protgroup2ac_df, proteins_df, iactions_df,
      protac_colls, pg_colls, terms_df, goterm_info_df,
      pcomplexes_df, pcomplex_iactors_df)

else

@load(joinpath(scratch_path, "$(proj_info.id)_omics_collections_$(proj_info.oesc_ver).jld2"),
      protgroups_df, protgroup2ac_df, proteins_df, iactions_df,
      protac_colls, pg_colls, terms_df, goterm_info_df,
      pcomplexes_df, pcomplex_iactors_df)

ProtgroupType = eltype(iactions_df.protgroup_id)

end

pg_id2name = sizehint!(Dict{ProtgroupType, String}(), nrow(protgroups_df));
for r in eachrow(protgroups_df)
    pg_id = r.protgroup_id
    if !ismissing(pg_id) && !ismissing(r.gene_names)
        pg_id2name[pg_id] = r.gene_names
    end
end

protgroup_iact_data_df = unique(iactions_df[!, [:protgroup_id, :single_uniprot_id]])
pg_id2ac = sizehint!(Dict{ProtgroupType, String}(), nrow(protgroup_iact_data_df));
for r in eachrow(protgroup_iact_data_df)
    pg_id = r.protgroup_id
    if !ismissing(pg_id) && !ismissing(r.single_uniprot_id)
        pg_id2ac[pg_id] = r.single_uniprot_id
    end
end

using OptEnrichedSetCover
include(joinpath(misc_scripts_path, "optcover_utils.jl"));

@info "Preparing mosaics..."
observed_protacs = Set(protgroup2ac_df[protgroup2ac_df.is_majority, :protein_ac])
pg_mosaics = OptCoverUtils.collections2mosaics(pg_colls,
                                      protac_colls, observed_protacs,
                                      setXset_frac_extra_elms=0.05,
                                      verbose=true);

BaitType = CategoricalArrays.leveltype(eltype(iactions_df[!, bait_col]))
pg_vorf_iactors = Dict{BaitType, Set{ProtgroupType}}()
for vorf_subdf in groupby(iactions_df, bait_col)
    pg_vorf_iactors[vorf_subdf[1, bait_col]] = Set(vorf_subdf.protgroup_id)
end

pg_iaction_mosaics = Dict(begin
    @info "Masking $mosaic_name by AP/MS interactions..."
    mosaic_name => OptCoverUtils.automask(mosaic, pg_vorf_iactors, verbose=true,
                                          max_sets=2000, min_nmasked=2, max_setsize=200, max_pvalue=1E-3)
end for (mosaic_name, mosaic) in pairs(pg_mosaics));

@save(joinpath(scratch_path, "$(proj_info.id)_iaction_mosaics_$(proj_info.oesc_ver).jld2"),
      proj_info,
      pg_mosaics, pg_iaction_mosaics)

using OptEnrichedSetCover

cover_params = CoverParams(setXset_factor=0.4,
                           uncovered_factor=0.0, covered_factor=0.0)
#cover_params = CoverParams(sel_prob=0.05, setXset_factor=0.5,
#                           uncovered_factor=0.5, covered_factor=0.01)

ENV["MKL_NUM_THREADS"] = 1

pg_iaction_mosaics_v = collect(pairs(pg_iaction_mosaics))
pg_iaction_covers_v = similar(pg_iaction_mosaics_v, Pair)
Threads.@threads for i in eachindex(pg_iaction_mosaics_v)
    mosaic_name, masked_mosaic = pg_iaction_mosaics_v[i]
    @info "Covering $mosaic_name by AP/MS interactions..."
    pg_iaction_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=2_000_000, WeightDigits=2,
                                    NWorkers=1,#Threads.nthreads()-1,
                                    MaxRestarts=200),
            true)
end
pg_iaction_covers = Dict(k => v for (k, v) in pg_iaction_covers_v)

using Printf

pg_iaction_covers_df = OptCoverUtils.covers_report(
    pg_iaction_covers, pg_vorf_iactors,
    pg_colls, pg_mosaics,
    pg_id2name, terms_df,
    maskid_col=[:bait_vorf],
    maskedset_col_prefix="bait");
pg_iaction_covers_df.bait_label =
    [isinteger(lbl) ? @sprintf("%.0f", lbl) : @sprintf("%.1f", lbl)
     for lbl in pg_iaction_covers_df.bait_vorf];
pg_iaction_covers_df.is_proteincomplex =
    (pg_iaction_covers_df.term_id .∈ Ref(Set(string.(proteincomplexes_goterms)))) .|
    (pg_iaction_covers_df.coll_id .== "protein_complexes");
pg_iaction_covers_df = filter!(r -> coalesce(r.intersect_genes, "") != "", pg_iaction_covers_df)

# reports that use protgroup_ids instead of gene names
pg_iaction_covers_pg_df = OptCoverUtils.covers_report(
    pg_iaction_covers, pg_vorf_iactors,
    pg_colls, pg_mosaics,
    nothing, terms_df,
    maskid_col=[:bait_vorf],
    maskedset_col_prefix="bait",
    objs_col_suffix="protgroups", maxels=nothing);
pg_iaction_covers_pg_df.is_proteincomplex =
    (pg_iaction_covers_pg_df.term_id .∈ Ref(Set(string.(proteincomplexes_goterms)))) .|
    (pg_iaction_covers_pg_df.coll_id .== "protein_complexes");
pg_iaction_covers_pg_df.bait_label =
    [isinteger(lbl) ? @sprintf("%.0f", lbl) : @sprintf("%.1f", lbl)
     for lbl in pg_iaction_covers_pg_df.bait_vorf];
pg_iaction_covers_pg_df = filter!(r -> r.intersect_protgroups != "", pg_iaction_covers_pg_df)
pg_iaction_covers_pg_long_df = DelimDataUtils.expand_delim_column(
    filter(r -> (r.nmasked >= 3) && (r.set_overlap_log10pvalue <= 1E-3), pg_iaction_covers_pg_df),
    list_col=:intersect_protgroups, elem_col=:protgroup_id,
    key_col=[:bait_label, :coll_id, :term_id, :term_name, :term_descr, :is_proteincomplex],
    delim=" ")
pg_iaction_covers_pg_long_df.protgroup_id = parse.(Int, pg_iaction_covers_pg_long_df.protgroup_id)
pg_iaction_covers_pg_long_df.protein_ac = get.(Ref(pg_id2ac), pg_iaction_covers_pg_long_df.protgroup_id, missing)
pg_iaction_covers_pg_data_df = combine(groupby(pg_iaction_covers_pg_long_df,
                                       [:coll_id, :protgroup_id, :protein_ac])) do pg_terms_df
    unique_terms_df = sort!(unique!(pg_terms_df[:, [:term_id, :term_name, :term_descr, :is_proteincomplex]]), :term_id)
    DataFrame(bait_labels = join(sort!(unique(pg_terms_df.bait_label)), '|'),
              term_ids = join(unique_terms_df.term_id, '|'),
              term_names = join(replace!(unique_terms_df.term_name, missing=>""), '|'),
              term_descrs = join(replace!(unique_terms_df.term_descr, missing=>""), '|'),
              is_proteincomplex = any(unique_terms_df.is_proteincomplex))
end

for pg_terms_df in groupby(pg_iaction_covers_pg_data_df, :coll_id)
    coll_id = pg_terms_df.coll_id[1]
    report_df = rename!(select!(pg_terms_df[:, :], Not(:coll_id)),
                        :term_ids => Symbol(coll_id, "_term_ids"),
                        :term_names => Symbol(coll_id, "_term_names"),
                        :term_descrs => Symbol(coll_id, "_term_descrs"))
    CSV.write(joinpath(analysis_path, "reports", "for_network", "$(proj_info.id)_pg_groups_oesc_$(proj_info.oesc_ver)_$(coll_id).txt"),
              report_df, missingstring="", delim='\t');
end

multicover_fname = "$(proj_info.id)_oesc_vorf_multicover_$(proj_info.oesc_ver).jld2"
@save(joinpath(scratch_path, multicover_fname),
      pg_colls, pg_iaction_mosaics, pg_iaction_covers, pg_iaction_covers_df)

if !@isdefined(pg_iaction_multicover_colls)
@load(joinpath(scratch_path, multicover_fname),
      pg_colls, pg_iaction_mosaics, pg_iaction_covers, pg_iaction_covers_df)
end;

pg_iaction_covers_signif_df = combine(groupby(pg_iaction_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df[1, :term_collection])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:bait_vorf, :bait_label],
                                                   max_term_pvalue=1E-3, max_set_pvalue=1E-2, min_set_overlap=nothing),
                   Not(:term_collection))
end

# expand all protein-protein interactions for enriched protein complexes (GO_CC + IntAct/CORUM)
pg_cocomplex_pairs_df = combine(groupby(pg_iaction_covers_pg_long_df[pg_iaction_covers_pg_long_df.is_proteincomplex, :],
                                        :coll_id)) do pg_coll_df
    coll_prepairs_df = combine(groupby(pg_coll_df, :term_id)) do pg_term_df
        pgs_df = pg_term_df[!, [:protgroup_id, :term_id]]
        pairs_df = innerjoin(rename(pgs_df, :protgroup_id => :src_protgroup_id),
                             rename(pgs_df, :protgroup_id => :dest_protgroup_id),
                             on=:term_id)
        filter!(r -> r.src_protgroup_id < r.dest_protgroup_id, pairs_df)
        return pairs_df
    end
    return combine(groupby(coll_prepairs_df, [:src_protgroup_id, :dest_protgroup_id])) do pair_df
        term_ids = Set(pair_df.term_id)
        unique_terms_df = filter(r -> r.term_id ∈ term_ids, terms_df)
        DataFrame(term_ids = join(unique_terms_df.term_id, '|'),
                  term_names = join(filter(!ismissing, unique_terms_df.term_name), '|'),
                  term_descrs = join(filter(!ismissing, unique_terms_df.term_descr), '|'))
    end
end
CSV.write(joinpath(analysis_path, "reports", "for_network", "$(proj_info.id)_pg_cocomplex_pairs_oesc_$(proj_info.oesc_ver).txt"),
          pg_cocomplex_pairs_df, missingstring="", delim='\t');

using PlotlyJS, TextWrap, Printf

include(joinpath(base_scripts_path, "misc_jl", "optcover_plots.jl"))
include(joinpath(base_scripts_path, "misc_jl", "optcover_heatmap.jl"))

oesc_plots_path = joinpath(plots_path, "oesc_$(proj_info.oesc_ver)")
isdir(oesc_plots_path) || mkdir(oesc_plots_path)
pareto_plots_path = joinpath(oesc_plots_path, "pareto")
isdir(pareto_plots_path) || mkdir(pareto_plots_path)
for (plot_mosaic, cover_coll) in pg_iaction_covers
    isempty(cover_coll.results) && continue
    paretofront_plot = OptCoverPlots.plot_paretofront(cover_coll.results[1], plot_unfolded=true)
    @info "Plotting $plot_mosaic Pareto front"
    #savefig(paretofront_plot,
    #        joinpath(plots_path, "$(proj_info.id)_$(plot_mosaic)_pareto.svg"))
    PlotlyJS.savehtml(paretofront_plot,
            joinpath(pareto_plots_path, "$(proj_info.id)_$(plot_mosaic)_oesc_$(proj_info.oesc_ver)_paretofront.html"))
end

function process_bait_axis(baits_df)
    bait_labels = baits_df.bait_label

    return baits_df, "v" .* bait_labels, bait_labels
end

heatmap_layout_attrs = Dict(
    ("SigDB_C2", true) => Dict(:margin_l => 500),
    ("SigDB_C2", false) => Dict(:margin_l => 500),
    ("Reactome", true) => Dict(:margin_l => 500),
    ("Reactome", false) => Dict(:margin_l => 500),
    ("Reactome", true) => Dict(:margin_l => 600),
    ("Reactome", false) => Dict(:margin_l => 600),
    ("GOcc", true) => Dict(:margin_l => 200),
    ("GOcc", false) => Dict(:margin_l => 200),
)

for term_coll in unique(pg_iaction_covers_df.term_collection), signif in (false, true)
    @info "Plotting $(signif ? "signif " : "")bait vORF heatmap for $term_coll..."
    layout_attrs = get(heatmap_layout_attrs, (term_coll, signif), Dict())
    df = filter(r -> r.term_collection == term_coll, signif ? pg_iaction_covers_signif_df : pg_iaction_covers_df)
    if nrow(df) == 0
        @warn "No term_collection=$term_coll rows"
        continue
    end

    coll_heatmap = OptCoverHeatmap.oesc_heatmap(df,
            elements_label="protein",
            mask_axis_title = "vORF",
            mask_cols = [:bait_label, :nbait],
            process_mask_axis=process_bait_axis,
            process_term_axis=OptCoverHeatmap.process_term_axis,
            margin_l=get(layout_attrs, :margin_l, 400),
            margin_b=get(layout_attrs, :margin_b, 80),
            transpose=false,
            colorscale = "Hot", reversescale=false,
            plot_bgcolor="#FFF", gridcolor="#DDD",#outformat in ["svg", "pdf"] ? "#000" : "#BBB",
            cell_width=25, cell_height=20, gridwidth=1)
    if coll_heatmap === nothing
        @warn "Empty heatmap for $term_coll"
    end
    for (k, v) in [#:width=>800, :height=>400,
                   :margin_r=>80,
                   :yaxis_tickfont_size=>12, :xaxis_tickangle=>45]
        coll_heatmap.plot.layout[k] = v
    end
    for outformat in ["svg", "pdf", "html"]
        plot_fname = joinpath(oesc_plots_path,
            "$(proj_info.id)_$(proj_info.oesc_ver)_$(term_coll)_X_vorfs$(signif ? "_signif" : "")_heatmap.$(outformat)")
        try
            if outformat == "html"
                PlotlyJS.savehtml(coll_heatmap, plot_fname, :embed);
            else
                savefig(coll_heatmap, plot_fname, width=coll_heatmap.plot.layout[:width], height=coll_heatmap.plot.layout[:height]);
            end
        catch e
            @warn "$term_coll generation failed: $e"
        end
    end
end

xlsx_compat(col::AbstractVector) =
    nonmissingtype(eltype(col)) <: Symbol || nonmissingtype(eltype(col)) <: CategoricalValue ?
    ifelse.(ismissing.(col), col, String.(col)) : col

using XLSX
include(joinpath(misc_scripts_path, "xlsx_utils.jl"));

XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_pg_vorfs_oesc_$(proj_info.oesc_ver).xlsx"),
                XLSXUtils.xlsx_compat(filter(r -> r.nmasked > 0, pg_iaction_covers_df)), overwrite=true, sheetname="enrichment report")
XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_pg_vorfs_oesc_signif_$(proj_info.oesc_ver).xlsx"),
                XLSXUtils.xlsx_compat(filter(r -> r.nmasked > 0, pg_iaction_covers_signif_df)), overwrite=true, sheetname="enrichment report")

@info "Covering the set of all interactors taken together..."
pg_all_iactors = Dict("all" => Set(iactions_df.protgroup_id))

pg_all_iactors_mosaics = Dict(begin
    @info "Masking $mosaic_name by all AP/MS interactors..."
    mosaic_name => OptCoverUtils.automask(mosaic, pg_all_iactors, verbose=true,
                                          max_sets=2000, min_nmasked=2, max_setsize=2000, max_pvalue=1E-3)
end for (mosaic_name, mosaic) in pairs(pg_mosaics));

pg_all_iactors_mosaics_v = collect(pairs(pg_all_iactors_mosaics))
pg_all_iactors_covers_v = similar(pg_all_iactors_mosaics_v, Pair)
Threads.@threads for i in eachindex(pg_all_iactors_mosaics_v)
    mosaic_name, masked_mosaic = pg_all_iactors_mosaics_v[i]
    @info "Covering $mosaic_name by all AP/MS interactors..."
    pg_all_iactors_covers_v[i] =
        mosaic_name => collect(masked_mosaic, cover_params,
            CoverEnumerationParams(max_set_score=0.0, max_covers=1),
            MultiobjOptimizerParams(ϵ=[0.02, 0.2], MaxSteps=5_000_000, WeightDigits=2,
                                    NWorkers=1,#Threads.nthreads()-1,
                                    MaxRestarts=200),
            true)
end
pg_all_iactors_covers = Dict(k => v for (k, v) in pg_all_iactors_covers_v)

#alt_cover_params = CoverParams(setXset_factor=0.5,
#                               uncovered_factor=0.0, covered_factor=0.0)

pg_all_iactors_covers_df = OptCoverUtils.covers_report(
    pg_all_iactors_covers, pg_all_iactors,
    pg_colls, pg_mosaics,
    pg_id2name, terms_df,
    #cover_params=alt_cover_params,
    maskid_col=[:dummy],
    maskedset_col_prefix="iactors")

pg_all_iactors_covers_pg_df = OptCoverUtils.covers_report(
    pg_all_iactors_covers, pg_all_iactors,
    pg_colls, pg_mosaics,
    nothing, terms_df,
    maskid_col=[:dummy],
    maskedset_col_prefix="iactors",
    objs_col_suffix="protgroups", maxels=nothing);
pg_all_iactors_covers_pg_df = select!(filter!(r -> r.intersect_protgroups != "", pg_all_iactors_covers_pg_df),
                                      Not(:dummy)) # there's only single mask

pg_all_iactors_covers_pg_long_df = DelimDataUtils.expand_delim_column(
    filter(r -> (r.nmasked >= 3) && (r.set_overlap_log10pvalue <= 1E-3), pg_all_iactors_covers_pg_df),
    list_col=:intersect_protgroups, elem_col=:protgroup_id,
    key_col=[:coll_id, :term_id, :term_name, :term_descr, :set_overlap_log10pvalue],
    delim=" ")
pg_all_iactors_covers_pg_long_df.protgroup_id = parse.(Int, pg_all_iactors_covers_pg_long_df.protgroup_id)
pg_all_iactors_covers_pg_long_df.protein_ac = get.(Ref(pg_id2ac), pg_all_iactors_covers_pg_long_df.protgroup_id, missing)

pg_all_iactors_covers_terms_df = combine(groupby(pg_all_iactors_covers_pg_long_df,
                                           [:coll_id, :term_id, :term_name, :term_descr])) do term_df
    DataFrame(protein_acs = join(term_df.protein_ac, ' '),
              protgroup_ids = join(term_df.protgroup_id, ' '),
              n_protgroups = nrow(term_df),
              set_overlap_log10pvalue = term_df.set_overlap_log10pvalue[1])
end

XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_pg_alliactors_protgroups_$(proj_info.oesc_ver).xlsx"),
                collect(xlsx_compat.(eachcol(pg_all_iactors_covers_pg_long_df))),
                DataFrames.names(pg_all_iactors_covers_pg_long_df), overwrite=true);
XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_pg_alliactors_terms_$(proj_info.oesc_ver).xlsx"),
                collect(xlsx_compat.(eachcol(pg_all_iactors_covers_terms_df))),
                DataFrames.names(pg_all_iactors_covers_terms_df), overwrite=true);

allcover_fname = "$(proj_info.id)_oesc_alliactors_cover_$(proj_info.oesc_ver).jld2"
@save(joinpath(scratch_path, allcover_fname),
      pg_colls, pg_all_iactors_mosaics, pg_all_iactors_covers, pg_all_iactors_covers_df)

if !@isdefined(pg_all_iactors_cover_colls)
@load(joinpath(scratch_path, allcover_fname),
      pg_colls, pg_all_iactors_mosaics, pg_all_iactors_covers, pg_all_iactors_covers_df)
end;

pg_all_iactors_covers_signif_df = combine(groupby(pg_all_iactors_covers_df, :term_collection)) do coll_df
    @info "Processing $(coll_df[1, :term_collection])..."
    return select!(OptCoverUtils.filter_multicover(coll_df, set_cols=[:dummy],
                                                   max_term_pvalue=1E-3, max_set_pvalue=1E-2, min_set_overlap=nothing),
                   Not(:term_collection))
end
select!(pg_all_iactors_covers_df, Not(:dummy)); # there's only single mask
select!(pg_all_iactors_covers_signif_df, Not(:dummy)); # there's only single mask

XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_pg_alliactors_covers_$(proj_info.oesc_ver).xlsx"),
                collect(xlsx_compat.(eachcol(filter(r -> r.nmasked > 0, pg_all_iactors_covers_df)))),
                DataFrames.names(pg_all_iactors_covers_df), overwrite=true);
XLSX.writetable(joinpath(analysis_path, "reports", "$(proj_info.id)_pg_alliactors_covers_signif_$(proj_info.oesc_ver).xlsx"),
                collect(xlsx_compat.(eachcol(filter(r -> r.nmasked > 0, pg_all_iactors_covers_signif_df)))),
                DataFrames.names(pg_all_iactors_covers_signif_df), overwrite=true);

sel_protgroups_df = protgroups_df[occursin.(Ref(r"RALGAP"), coalesce.(protgroups_df.gene_names, "")), :]
sel_protgroups_df = protgroups_df[occursin.(Ref(r"MOV10|USP11|NKIRAS2"), coalesce.(protgroups_df.gene_names, "")), :]
sel_protgroup_ids = sel_protgroups_df.protgroup_id

pg_iaction_mosaics[:GObp]
sel_masks_ixs = findall(keys(pg_vorf_iactors) .∈ Ref(Set(["66", "66c"])))

sel_term_ids = Dict(coll_id => collect(keys(coll))[[any(pg_id -> pg_id ∈ sel_protgroup_ids, pg_ids) for (term_id, pg_ids) in coll]]
    for (coll_id, coll) in pg_colls);
sel_mosaic_ids = Dict(begin
    pg_mosaic = pg_iaction_mosaics[coll_id]
    setixs = [pg_mosaic.original.set2ix[term_id] for term_id in term_ids]
    filter!(setix -> begin
        haskey(pg_mosaic.set2masks, setix) || return false
        any(molap -> molap.mask ∈ sel_masks_ixs, pg_mosaic.set2masks[setix])
    end, setixs)
    coll_id => setixs
end for (coll_id, term_ids) in sel_term_ids)

sel_problem = MultiobjCoverProblem(pg_iaction_mosaics[:GObp])
sel_vars = findall(setix -> setix ∈ sel_mosaic_ids[:GObp], sel_problem.var2set)
sum(sel_problem.nmasked_tile[sel_problem.tileXvar[:, 122]])
sum(sel_problem.nmasked_tile[sel_problem.tileXvar[:, 1017]])
sel_weights = pg_iaction_multicover_colls[:GObp].results[1].weights
@show OptEnrichedSetCover.score(sel_problem, sel_weights)
sel_vars = [818]# findall(w -> w > 0, sel_weights)
for varix in sel_vars
    varw = fill(0.0, length(sel_weights))
    varw[varix] = 1.0
    setXset = fill!(similar(varw), 0.0)
    OptEnrichedSetCover.minplus_bilinear!(min, setXset, sel_problem.varXvar_scores, sel_weights, varw)
    @show varix OptEnrichedSetCover.miscover_score(sel_problem, varw) OptEnrichedSetCover.score(sel_problem, varw)
    @show sum(setXset)
end
sel_problem.tileXvar[:, 818]
sel_problem.var2set[818]
pg_iaction_mosaics[:GObp].original.ix2set[4780]
sel_problem.varXvar_scores[818, :]

OptEnrichedSetCover.miscover_score(sel_problem, pg_iaction_multicover_colls[:GObp].results[1].weights)

sel_set_info_df = omics_set_info_df[[], :]
for (coll_id, pg_coll) in pg_colls
    @info "Searching $coll_id"
    for (set_id, pg_set) in pg_coll
        if any(x -> in(x, pg_set), sel_protgroup_ids)
            set_info_df = omics_set_info_df[(omics_set_info_df.coll_id .== string(coll_id)) .&
                                            (omics_set_info_df.term_id .== set_id), :]
            if isempty(set_info_df)
                @warn set_id
            else
                append!(sel_set_info_df, set_info_df)
            end
        end
    end
end
sel_set_info_df

OptCoverUtils.automask(pg_mosaics[:GObp], values(pg_vorf_iactors), max_sets=1500, min_nmasked=2)
