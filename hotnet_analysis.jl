proj_info = (id = "vgirault_vzvapms",
             data_ver = "20180301",
             analysis_ver = "20210615",
             hotnet_ver = "20210707",
             model_obj = "protgroup")
using Pkg
Pkg.activate(joinpath(base_scripts_path, "adhoc", proj_info.id))

using Revise
using Distances, DataFrames, DataFramesMeta, CSV, RData, CodecZstd, JLD2

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl")
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(analysis_path, "scratch")
const plots_path = joinpath(analysis_path, "plots")

includet(joinpath(misc_scripts_path, "frame_utils.jl"))
includet(joinpath(misc_scripts_path, "delimdata_utils.jl"))
includet(joinpath(misc_scripts_path, "graphml_writer.jl"))
includet(joinpath(misc_scripts_path, "forceatlas3_layout.jl"))
includet(joinpath(misc_scripts_path, "hotnet_utils.jl"))
includet(joinpath(misc_scripts_path, "fasta_reader.jl"))

using Revise, DataFrames, CategoricalArrays, CSV, RData, SimpleWeightedGraphs, Graphs,
      Statistics, StatsBase, LinearAlgebra, Chain
using HierarchicalHotNet
HHN = HierarchicalHotNet

# Import reactome
reactomefi_diedges_df = HotnetUtils.import_reactomefi(joinpath(party3rd_data_path, "FIsInGene_122220_with_annotations.txt"), verbose=true)
reactomefi_genes = levels(reactomefi_diedges_df.gene1)
reactomefi_gene2index = Dict(gene => ix for (ix, gene) in enumerate(reactomefi_genes))

# reverse-direction graph
reactomefi_digraph_revfull, reactomefi_digraph_revfull_vertices =
    HHN.import_digraph(reactomefi_diedges_df, src_col=:gene2, dest_col=:gene1, weight_col=:score)
@assert levels(reactomefi_diedges_df.gene1) == reactomefi_digraph_revfull_vertices
# find non-tirivial connected (not necessarily strongly) components
reactomefi_digraph_revfull_conncomp = HHN.strongly_connected_components(
    Graphs.adjacency_matrix(reactomefi_digraph_revfull) + transpose(Graphs.adjacency_matrix(reactomefi_digraph_revfull)),
    HHN.EdgeTest{Float64}(threshold=nothing))
reactomefi_digraph_revfull_conncomp_used = filter(comp -> length(comp) > 1, reactomefi_digraph_revfull_conncomp)

reactomefi_genes_mask = eachindex(reactomefi_digraph_revfull_vertices) .∈ Ref(reactomefi_digraph_revfull_conncomp_used.elems)
@info "$(sum(reactomefi_genes_mask)) gene(s) used"
vertex2gene = findall(reactomefi_genes_mask)
gene2vertex = fill(0, length(reactomefi_genes))
@inbounds for (vertex, gene) in enumerate(vertex2gene)
    gene2vertex[gene] = vertex
end

# reversed reactome FI subgraph induced by reactomefi_genes_mask
reactomefi_diedges_df[!, :vertex1] = gene2vertex[levelcode.(reactomefi_diedges_df.gene1)]
reactomefi_diedges_df[!, :vertex2] = gene2vertex[levelcode.(reactomefi_diedges_df.gene2)]
reactomefi_diedges_used_df = filter(r -> (r.vertex1 > 0) && (r.vertex2 > 0), reactomefi_diedges_df)
reactomefi_digraph_rev, reactomefi_digraph_vertex_indices =
    HHN.import_digraph(reactomefi_diedges_used_df,
                       src_col=:vertex2, dest_col=:vertex1, weight_col=:score)
@assert reactomefi_digraph_vertex_indices == 1:nv(reactomefi_digraph_rev)


# load experiments data
analysis_rdata = load(joinpath(analysis_path, "results", "$(proj_info.id)_apmsXoefp_$(proj_info.analysis_ver).RData"))

#--- map OE proteome to reactome network
objects_df = copy(analysis_rdata["oefp_protgroups.df"])
objects_df.is_used = .!objects_df.is_contaminant .& .!objects_df.is_viral
oeproteome_contrasts_df = @chain analysis_rdata["protgroupXbait_tests.df"] begin
    semijoin(_, filter(r -> r.is_used, objects_df), on=:protgroup_id)
    @rtransform!(_, :p_value = :p_value_wilcox, :p_value_adj = :p_value_wil_BH)
    rename!(_, :protgroup_id => :oefp_protgroup_id)
    @aside countmap(_.condition)
    # calculate vertex weights
    @rtransform!(_,
        :is_source =
            (coalesce(:p_value_adj, 1.0) <= 0.05 && abs(:log2_foldchange) >= 1.0),
        :vertex_weight = ifelse(:is_source,
            (-log10(max(1E-20, :p_value)))^0.5 * abs(:log2_foldchange)^0.5, 0.0),
        :edge_mult = 1.0 + ifelse(
                (coalesce(:p_value, 1.0) <= 0.001 && abs(:log2_foldchange) >= 0.25),
                (-log10(max(1E-20, :p_value)))^0.5 * abs(:log2_foldchange)^0.5, 0.0))
    @aside extrema(_[_.vertex_weight .> 0, :vertex_weight])
    @aside quantile(_[_.vertex_weight .> 0, :vertex_weight], 0.5)
end

oeproteome_contrasts_stats_df = combine(groupby(oeproteome_contrasts_df, [:msbatch, :condition, :contrast]),
        :vertex_weight => (w -> any(>(0), w) ? quantile(filter(>(0), w), 0.5) : missing) => :vertex_weight_median,
        :vertex_weight => (w -> sum(>(0), w)) => :n_nzweights,
        :is_source => (w -> sum(w)) => :n_sources)
flows_path = joinpath(analysis_path, "networks", "apms_oeproteome_flows_$(proj_info.hotnet_ver)")
isdir(flows_path) || mkdir(flows_path)
CSV.write(joinpath(flows_path, "oeproteome_contrast_stats_$(proj_info.hotnet_ver).txt"),
          oeproteome_contrasts_stats_df, delim='\t')

proteins_df = vcat(begin
    host_df = Fasta.read_uniprot(joinpath(analysis_path, "fasta", "uniprot-9606_proteome_human_reviewed_canonical_isoforms_191008.fasta"))
    host_df[!, :is_viral] .= false
    host_df
end, begin
    viral_df = Fasta.read_uniprot(joinpath(analysis_path, "fasta", "VZV_vOka_20171208.fasta"))
    viral_df[!, :is_viral] .= true
    viral_df
end)
proteins_df[!, :is_contaminant] .= false

oeproteome_obj2gene_df = innerjoin(
    filter(r -> r.is_majority, analysis_rdata["oefp_protgroup2protein.df"]),
    filter(r -> !ismissing(r.genename), proteins_df),
    on=:protein_ac)
oeproteome_obj2gene_df[!, :gene_id] = get.(Ref(reactomefi_gene2index), oeproteome_obj2gene_df.genename, 0)
oeproteome_gene_effects_df = combine(groupby(innerjoin(oeproteome_contrasts_df, oeproteome_obj2gene_df, on=:oefp_protgroup_id),
          [:condition, :genename, :gene_id])) do df
    rowix = findmin(df.p_value)[2]
    return df[rowix:rowix, :]
end;
countmap(oeproteome_gene_effects_df.vertex_weight .> 0)

geneinfo_df = combine(groupby(filter(r -> !ismissing(r.genename), proteins_df), [:genename])) do gene_df
                  gene_df[1, :]
              end

bait2weights = Dict{String, NamedTuple{(:vertex, :edge), Tuple{Vector{Float64}, Vector{Float64}}}}()
for bait_genes_df in groupby(filter(r -> r.gene_id > 0,
                                    oeproteome_gene_effects_df), :condition)
    vertex_weights = fill(0.0, nv(reactomefi_digraph_rev))
    edge_mults = fill(1.0, nv(reactomefi_digraph_rev))
    for r in eachrow(bait_genes_df)
        v = gene2vertex[r.gene_id]
        v2 = v > 0 ? searchsortedfirst(reactomefi_digraph_vertex_indices, v) : 0
        @assert v2 == v "v2=$v2 doesn't match v=$v"
        if v > 0
            vertex_weights[v] = r.vertex_weight
            edge_mults[v] = r.edge_mult
        end
    end
    bait2weights[bait_genes_df.condition[1]] = (vertex=vertex_weights, edge=edge_mults)
end
bait_ids = sort!(collect(keys(bait2weights)))

#region map APMS data to reactome network
apms_iactions_df = copy(analysis_rdata["vzv_apmsXoefp.df"])
# remove upregulated proteins that are potential false positive APMS interactions
filter!(r -> !coalesce(r.oefp_is_upreg, false), apms_iactions_df)

apms_baits_df = unique!(select(apms_iactions_df, [:display_bait, :bait_id]))
oefp_conditions_df = unique!(select(oeproteome_contrasts_df, :condition))

apmsXoefp_conditions_df = CSV.read(joinpath(analysis_path, "networks", "apmsXoefp_pairs_$(proj_info.analysis_ver).txt"), DataFrame, delim='\t')

flows_path = joinpath(analysis_path, "networks", "apms_oeproteome_flows_$(proj_info.hotnet_ver)")
isdir(flows_path) || mkdir(flows_path)

apms_obj2gene_df = innerjoin(
    filter(r -> r.is_majority, analysis_rdata["apms_protgroup2protein.df"]),
    filter(r -> !ismissing(r.genename), proteins_df),
    on=:protein_ac)
apms_obj2gene_df.apms_protgroup_id = convert(Vector{Int}, apms_obj2gene_df.apms_protgroup_id)
apms_obj2gene_df[!, :gene_id] = get.(Ref(reactomefi_gene2index), apms_obj2gene_df.genename, 0)
apms_gene_iactions_df = combine(groupby(innerjoin(apms_iactions_df, apms_obj2gene_df, on=:apms_protgroup_id),
          [:display_bait, :genename, :gene_id])) do df
    return df[findmin(df.p_value)[2], :]
end;
apms_gene_iactions_df = innerjoin(apms_gene_iactions_df,
        filter(r -> !ismissing(r.condition) && !ismissing(r.display_bait), apmsXoefp_conditions_df),
        on=[:bait_id, :display_bait])
apms_gene_iactions_df.is_sink = .!apms_gene_iactions_df.oefp_is_upreg

bait2apms_genes = Dict{String, Vector{Int}}()
bait2apms_vertices = Dict{String, Vector{Int}}()
for bait_genes_df in groupby(filter(r -> r.is_sink, apms_gene_iactions_df), :condition)
    baitkey = bait_genes_df.condition[1]
    gene_ids = sort!(unique!(collect(filter(>(0), bait_genes_df.gene_id))))
    bait2apms_genes[baitkey] = gene_ids
    vertices = unique!(filter!(!=(0), gene2vertex[gene_ids]))
    v2s = map(vertices) do v
        v2 = searchsortedfirst(reactomefi_digraph_vertex_indices, v)
        @assert v2 == v
        return (v2 <= length(reactomefi_digraph_vertex_indices)) &&
            (reactomefi_digraph_vertex_indices[v2] == v) ? v2 : 0
    end
    filter!(>(0), v2s)
    bait2apms_vertices[baitkey] = v2s
end
#endregion

randomwalk_params = (restart_prob = 0.4,
                     source_weight_min = 1E-3, # any non-zero weight is a source
                     step_weight_min = 1E-6, # flow edge with the weight below that is reset to zero
                     walk_weight_min = 1E-5,
                     flow_weight_max = nothing,
                     flow_pvalue_max = 0.05) # flow edge with the p-value below that is filtered out
reactomefi_adjmtx = Matrix(Graphs.adjacency_matrix(reactomefi_digraph_rev, dir=:in));
reactomefi_walkmtx = HHN.random_walk_matrix(reactomefi_digraph_rev, randomwalk_params.restart_prob)

#region select optimal restart probability
#=
test_restart_probs = 0.05:0.05:0.95
test_diedge_weight_min = 1.0
ref_bait = "VZV-12"
reactomefi_test_walkmtxs = [HHN.similarity_matrix(
        HHN.stepmatrix(reactomefi_adjmtx, inedge_weights=bait2weights[ref_bait].edge),
        bait2weights[ref_bait].vertex, restart_probability=p)
    for p in test_restart_probs]
=#
#endregion

using SparseArrays

sparseclamp(mtx::SparseMatrixCSC; absminval::Union{Number, Nothing} = nothing) = mtx
function sparseclamp(mtx::Matrix; absminval::Union{Number, Nothing} = nothing)
    if !isnothing(absminval)
        @inbounds for i in eachindex(mtx)
            if abs(mtx[i]) < absminval
                mtx[i] = 0
            end
        end
    end
    return sparse(mtx)
end

sel_bait_ids = bait_ids
bait2mtxs = Dict(bait_id => begin
    @info "Processing $bait_id"
    stepmtx = HHN.stepmatrix(reactomefi_adjmtx, inedge_weights=bait2weights[bait_id].edge)
    walkmtx = HHN.random_walk_matrix(stepmtx, randomwalk_params.restart_prob)
    vertex_weights = bait2weights[bait_id].vertex
    diff_vertex_weights = walkmtx * vertex_weights
    dwstepmtx = HHN.stabilized_stepmatrix(stepmtx, vertex_weights, randomwalk_params.restart_prob, diff_vertex_weights)
    dwstepmtx = sparseclamp(dwstepmtx, absminval=randomwalk_params.step_weight_min)

    (stepmtx = sparseclamp(stepmtx),
     walkmtx = sparseclamp(walkmtx, absminval=randomwalk_params.walk_weight_min),
     dwstepmtx = dwstepmtx,
     diedge_ixs = findall(!=(0), vec(dwstepmtx)),
     vertex_weights_diffused = diff_vertex_weights)
end for bait_id in sel_bait_ids)

# filter out zero matrices
bait2mtxs = Dict(bait_id => mtxs for (bait_id, mtxs) in pairs(bait2mtxs)
                 if sum(mtxs.vertex_weights_diffused) > 0 && nnz(mtxs.dwstepmtx) > 0)
valid_bait_ids = sort!(collect(keys(bait2mtxs)))

sum(>=(randomwalk_params.flow_weight_min), bait2mtxs["VZV-12"].dwstepmtx |> vec)
map(mtxs -> count(>=(randomwalk_params.flow_weight_min), nonzeros(mtxs.dwstepmtx)), values(bait2mtxs))
quantile(nonzeros(bait2mtxs["VZV-12"].dwstepmtx), 0.1)

sel_bait_ids = valid_bait_ids
reactomefi_trees = Vector{HHN.SCCTree{Float64}}(undef, length(sel_bait_ids));
Threads.@threads for i in eachindex(sel_bait_ids)
    bait_id = sel_bait_ids[i]
    walkmtx = bait2mtxs[bait_id].dwstepmtx
    @info "Bait $bait_id: $(size(walkmtx)) matrix"
    reactomefi_trees[i] = HHN.scctree(walkmtx, verbose=false, method=:bisect)
end
bait2tree = Dict(bait_id => reactomefi_trees[i] for (i, bait_id) in enumerate(sel_bait_ids))

# bin vertices mapped to pg_ids (i.e. quantified) according to their in/out degree
vertex_bins = HHN.vertexbins(reactomefi_digraph_rev, findall(>(-1), vertex2gene),
                             by=:outXin, method=:tree, nbins=100)
W = Float16 # permuted weight
bait2perm_weights = Dict{String, NamedTuple{(:vertex, :edge), Tuple{Matrix{W}, Matrix{W}}}}()
for (bait_id, weights) in bait2weights
    (bait_id ∈ valid_bait_ids) || continue
    vertex_weights_perms = fill(zero(W), length(weights.vertex), 1000)
    edge_weights_perms = fill(zero(W), length(weights.edge), 1000)
    vertex_perm = collect(eachindex(weights.vertex))
    for (vperm, eperm) in zip(eachcol(vertex_weights_perms), eachcol(edge_weights_perms))
        HHN.randpermgroups!(vertex_perm, vertex_bins)
        @inbounds for (i, j) in enumerate(vertex_perm)
            vperm[i] = weights.vertex[j]
            eperm[i] = weights.edge[j]
        end
    end
    bait2perm_weights[bait_id] = (vertex=vertex_weights_perms, edge=edge_weights_perms)
end

using Serialization, CodecZstd
open(ZstdCompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_prepare_$(proj_info.hotnet_ver).jlser.zst"), "w") do io
    serialize(io, (reactomefi_diedges_df, reactomefi_diedges_used_df,
              geneinfo_df, vertex2gene, gene2vertex, reactomefi_genes, reactomefi_digraph_rev,
              bait2apms_vertices, apms_gene_iactions_df, oeproteome_gene_effects_df,
              valid_bait_ids, bait2weights,
              randomwalk_params, bait2mtxs, bait2tree,
              vertex_bins, bait2perm_weights))
end

# save per-bait for use on a cluster
permtrees_input_prefix = "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver)"
sel_bait_ids = valid_bait_ids
open(ZstdCompressorStream, joinpath(scratch_path, "$(permtrees_input_prefix).jlser.zst"), "w") do io
    serialize(io, (vertex2gene, reactomefi_genes, reactomefi_digraph_rev,
                   [bait_id => size(bait2perm_weights[bait_id].vertex, 2) for bait_id in sel_bait_ids],
                   randomwalk_params))
end;

isdir(joinpath(scratch_path, permtrees_input_prefix)) || mkdir(joinpath(scratch_path, permtrees_input_prefix))
for (bait_ix, bait_id) in enumerate(sel_bait_ids)
    open(ZstdCompressorStream, joinpath(scratch_path, permtrees_input_prefix, "bait_$(bait_ix)_perms.jlser.zst"), "w") do io
        serialize(io, (bait_id, bait2weights[bait_id],
                       bait2perm_weights[bait_id],
                       get(bait2apms_vertices, bait_id, Int[]),
                       bait2mtxs[bait_id].stepmtx, bait2mtxs[bait_id].dwstepmtx, bait2mtxs[bait_id].walkmtx,
                       bait2mtxs[bait_id].diedge_ixs, bait2mtxs[bait_id].vertex_weights_diffused,
                       bait2tree[bait_id]))
    end
end;

using Serialization, CodecZstd
reactomefi_diedges_df, reactomefi_diedges_used_df,
geneinfo_df, vertex2gene, gene2vertex, reactomefi_genes,
reactomefi_digraph_rev,
bait2apms_vertices, apms_gene_iactions_df, oeproteome_gene_effects_df,
valid_bait_ids, bait2weights,
randomwalk_params, bait2mtxs, bait2tree,
vertex_bins, _ = open(deserialize, ZstdDecompressorStream,
                      joinpath(scratch_path, "$(proj_info.id)_hotnet_prepare_$(proj_info.hotnet_ver).jlser.zst"), "r")

using LinearAlgebra, HierarchicalHotNet
HHN = HierarchicalHotNet

# collect chunks of network diffusion / tree statistics
using Base.Filesystem, Serialization, CodecZstd

# collect real data-based statistics (see hotnet_treestats_chunk.jl)
chunk_prefix = "$(proj_info.id)_hotnet_treestats_$(proj_info.hotnet_ver)"
tree_stats_dfs = Vector{DataFrame}()
for (root, dirs, files) in Filesystem.walkdir(joinpath(scratch_path, chunk_prefix))
    if isempty(files)
        @warn "Found no files in $root"
        continue # skip the empty folder
    else
        @info "Found $(length(files)) file(s) in $root, processing..."
    end
    for file in files
        (startswith(file, chunk_prefix) && endswith(file, ".jlser.zst")) || continue
        chunk_proj_info, chunk_bait_ids, chunk_tree_stats_df =
            open(ZstdDecompressorStream, joinpath(root, file), "r") do io
            deserialize(io)
        end
        push!(tree_stats_dfs, chunk_tree_stats_df)
    end
end
tree_stats_df = vcat(tree_stats_dfs...)
tree_stats_dfs = nothing
@save(joinpath(scratch_path, "$(proj_info.id)_hotnet_treestats_$(proj_info.hotnet_ver).jld2"),
      tree_stats_df)

# collect assembled permuted data-based statistics (see hotnet_perm_chunk.jl and hotnet_perm_assemble.jl)
chunk_prefix = "$(proj_info.id)_hotnet_perm_assembled_$(proj_info.hotnet_ver)"
perm_tree_stats_dfs = Vector{DataFrame}()
vertex_stats_dfs = similar(perm_tree_stats_dfs)
diedge_stats_dfs = similar(perm_tree_stats_dfs)
chunks_update = Threads.SpinLock()
for (root, dirs, files) in Filesystem.walkdir(joinpath(scratch_path, chunk_prefix))
    if isempty(files)
        @warn "Found no files in $root"
        continue # skip the empty folder
    else
        @info "Found $(length(files)) file(s) in $root, processing..."
    end
    Threads.@threads for i in eachindex(files)
        file = files[i]
        (startswith(file, chunk_prefix) && endswith(file, ".jlser.zst")) || continue
        chunk_tree_stats_df, chunk_vertex_stats_df, chunk_diedge_stats_df =
            open(ZstdDecompressorStream, joinpath(root, file), "r") do io
            deserialize(io)
        end
        lock(chunks_update)
        push!(perm_tree_stats_dfs, chunk_tree_stats_df)
        push!(vertex_stats_dfs, chunk_vertex_stats_df)
        push!(diedge_stats_dfs, chunk_diedge_stats_df)
        unlock(chunks_update)
    end
end
perm_tree_stats_df = vcat(perm_tree_stats_dfs...)
vertex_stats_df = vcat(vertex_stats_dfs...)
diedge_stats_df = vcat(diedge_stats_dfs...)
perm_tree_stats_dfs = nothing
vertex_stats_dfs = nothing
diedge_stats_dfs = nothing
GC.gc()

vertex_stats_df = leftjoin(vertex_stats_df,
                           select(oeproteome_gene_effects_df, :condition => :bait_id, :gene_id, :is_source,
                                  :p_value => :oeproteome_p_value, :log2_foldchange => :oeproteome_median_log2),
                           on=[:bait_id, :gene_id])
vertex_stats_df2 = leftjoin(vertex_stats_df,
                            select(apms_gene_iactions_df, :condition => :bait_id, :genename => :gene_name,
                                   :p_value => :apms_p_value, :median_log2 => :apms_median_log2, :is_sink),
                            on=[:bait_id, :gene_name])
@assert nrow(vertex_stats_df) == nrow(vertex_stats_df2)
vertex_stats_df = vertex_stats_df2
vertex_stats_df.is_sink = coalesce.(vertex_stats_df.is_sink, false)
vertex_stats_df.is_hit = vertex_stats_df.prob_perm_walkweight_greater .<= 0.05
vertex_stats_df = leftjoin(vertex_stats_df, geneinfo_df, on=:gene_name=>:genename)

countmap(tree_stats_df.bait_id)
countmap(perm_tree_stats_df.bait_id)
@save(joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_stats_$(proj_info.hotnet_ver).jld2"),
      vertex_stats_df, diedge_stats_df, tree_stats_df, perm_tree_stats_df)

@load(joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_stats_$(proj_info.hotnet_ver).jld2"),
      vertex_stats_df, diedge_stats_df, tree_stats_df, perm_tree_stats_df)

tree_all_stats_df = perm_tree_stats_df
append!(tree_all_stats_df, select(tree_stats_df, Not(:threshold_bin)))
tree_all_stats_df.is_permuted = vcat(trues(nrow(tree_all_stats_df)-nrow(tree_stats_df)),
                                     falses(nrow(tree_stats_df)))
threshold_range = (0.0, 0.05)
tree_binstats_df = combine(groupby(tree_all_stats_df, :bait_id)) do bait_treestats_df
    @info "binstats($(bait_treestats_df.bait_id[1]))"
    if nrow(bait_treestats_df) == 0
        @warn "cutstats frame is empty, skipping"
        return DataFrame()
    end
    HHN.bin_treecut_stats(bait_treestats_df,
                          by_cols=[:is_permuted, :tree],
                          threshold_nbins=100, threshold_maxbinwidth=0.0001,
                          threshold_range=threshold_range)
end

# annotate tree_stats with threshold_bin indices
tree_stats_df2 = copy(tree_stats_df, copycols=false)
tree_stats_df2.threshold_bin = missings(Int, nrow(tree_stats_df2))
for df in groupby(tree_binstats_df, :bait_id)
    bait_bins = sort!([unique(collect(skipmissing(df.threshold_binmin))); [maximum(skipmissing(df.threshold_binmax))]])
    bait_stats_df = view(tree_stats_df2, tree_stats_df2.bait_id .== df.bait_id[1], :)
    HHN.add_bins!(bait_stats_df, :threshold, bait_bins)
end
tree_stats_df = tree_stats_df2

tree_perm_aggstats_df = combine(groupby(filter(r -> r.is_permuted, tree_binstats_df), :bait_id)) do perm_binstats_df
    @info "aggregate_treecut_binstats($(perm_binstats_df.bait_id[1]))"
    HHN.aggregate_treecut_binstats(perm_binstats_df, by_cols=[:is_permuted, :threshold_bin])
end
tree_perm_aggstats_wide_df = unstack(filter(r -> !ismissing(r.threshold_binmid), tree_perm_aggstats_df),
                                     [:bait_id, :is_permuted, :threshold_bin, :threshold_binmid], :quantile,
                                     intersect(HHN.TreecutMetrics, propertynames(tree_perm_aggstats_df)),
                                     namewidecols=(valcol, qtl, sep) -> Symbol(valcol, sep, HHN.quantile_suffix(qtl)))

using Distributions
nnodes_threshold_prior = Cauchy(200, 100)

cut_threshold_range = (1E-4, 1E-2)
tree_extremes_df = HHN.extreme_treecut_stats(
    filter(r -> !r.is_permuted && (cut_threshold_range[1] <= coalesce(r.threshold, -1.0) <= cut_threshold_range[2]), tree_stats_df),
    tree_perm_aggstats_df,
    threshold_weight = r -> pdf(nnodes_threshold_prior, r.topn_components_sizesum),
    relative = false,
    extra_join_cols = [:bait_id])

@save(joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_binstats_$(proj_info.hotnet_ver).jld2"),
      tree_stats_df, tree_binstats_df, tree_perm_aggstats_df, tree_perm_aggstats_wide_df, tree_extremes_df)

@load(joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_binstats_$(proj_info.hotnet_ver).jld2"),
      tree_stats_df, tree_binstats_df, tree_perm_aggstats_df, tree_perm_aggstats_wide_df, tree_extremes_df)

include(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))

sel_quantile = 0.25
bait2cut_threshold = Dict(begin
    bait_id = df[1, :bait_id]
    valtype = df[1, :value_type]
    thresh_df = nothing
    for (metric, stat) in [(:flow_avginvlen, "max"), (:topn_components_sizesum, "max"), (:maxcomponent_size, "max"),
                           (:flow_avgweight, "max"), (:flow_avghopweight, "max")]
        metric_df = filter(r -> (r.metric == metric) && (r.value_type == valtype) &&
                                (r.stat == stat) &&
                                (r.quantile == (stat == "max" ? 1 - sel_quantile : sel_quantile)), df)
        (nrow(metric_df) == 0) && continue
        @assert nrow(metric_df) == 1
        delta = metric_df[1, :delta]
        if !ismissing(delta) &&
           (((stat == "max") && (delta > 0.0)) ||
            ((stat == "min") && (delta < 0.0)))
            thresh_df = metric_df
            break
        end
    end
    if thresh_df !== nothing
        @info "$bait_id $valtype: using $(thresh_df.metric[1])=$(thresh_df.value[1]) (delta=$(thresh_df.delta[1])) for threshold=$(thresh_df.threshold[1])"
        thresh_df
    else
        @warn "No significant difference between real and permuted results for $bait_id"
    end
    (bait_id, valtype) => thresh_df
end for df in groupby(filter(r -> !ismissing(r.quantile), tree_extremes_df), [:bait_id, :value_type]))
#    if (r.type == "max") && !ismissing(r.flow_avgweight_perm_50))
bait_cut_thresholds_df = reduce(vcat, df for df in values(bait2cut_threshold) if !isnothing(df))
bait2treecut = Dict((baitid, valtype) =>
    !isnothing(threshold_df) ? HHN.cut(bait2tree[baitid], threshold_df[1, :threshold], minsize=1) : nothing
    for ((baitid, valtype), threshold_df) in bait2cut_threshold)

includet(joinpath(misc_scripts_path, "graphml_writer.jl"))

function flowgraphml(bait_id::AbstractString, threshold::Number; kwargs...)
    vertices_df = DataFrame(vertex = 1:nv(reactomefi_digraph_rev),
                            gene = ifelse.(vertex2gene .> 0, vertex2gene, missing),
                            gene_name = ifelse.(vertex2gene .> 0, getindex.(Ref(reactomefi_genes), vertex2gene), missing))
    bait_sinks = get(bait2apms_vertices, bait_id, Vector{Int}())
    vertices_df.apms_is_sink = vertices_df.vertex .∈ Ref(Set(bait_sinks))
    if @isdefined(vertex_stats_df)
        bait_vertex_stats_df = select!(filter(r -> r.bait_id == bait_id, vertex_stats_df),
            :vertex, :oeproteome_median_log2, :oeproteome_p_value, :is_source => :oeproteome_is_source,
            :apms_median_log2, :apms_p_value)
        vertices_df = leftjoin(vertices_df, bait_vertex_stats_df, on=:vertex)
        vertices_df.oeproteome_median_log2_signif = ifelse.(coalesce.(vertices_df.oeproteome_is_source, false),
                                                            vertices_df.oeproteome_median_log2, missing)
        vertices_df.apms_median_log2_signif = ifelse.(coalesce.(vertices_df.apms_is_sink, false),
                                                    vertices_df.apms_median_log2, missing)
    end

    # reverse original reactomefi diedges, since the flows are using the reversed ones
    orig_diedges = DataFrames.rename(reactomefi_diedges_used_df,
                            :vertex1 => :target, :vertex2 => :source,
                            :score => :weight, :direction => :interaction_type)
    HotnetUtils.fix_reactomefi_iactiontype!(orig_diedges.interaction_type)
    HotnetUtils.fix_reactomefi_iactiontype!(orig_diedges.diedge_type)
    baitwalk = bait2mtxs[bait_id]

    return HotnetUtils.flowgraphml(
            bait2tree[bait_id], threshold; step_threshold=threshold,
            walkmatrix=Matrix(baitwalk.dwstepmtx), stepmatrix=Matrix(baitwalk.dwstepmtx),
            vertices_labels=[v2g > 0 ? reactomefi_genes[v2g] : missing for v2g in vertex2gene],
            vertices_info=vertices_df,
            vertices_weights=bait2weights[bait_id].vertex,
            vertices_stats=@isdefined(vertex_stats_df) ? select!(filter(r -> r.bait_id == bait_id, vertex_stats_df),
                    Not([:gene_name, :oeproteome_median_log2, :oeproteome_p_value, :apms_median_log2, :apms_p_value, :is_source, :is_sink])) : nothing,
            diedges_stats=@isdefined(diedge_stats_df) ?
                select!(rename!(filter(r -> r.bait_id == bait_id, diedge_stats_df),
                        :src=>:source, :dest=>:target), Not([:weight, :walkweight]))
                : nothing,
            source_threshold=randomwalk_params.source_weight_min,
            edge_pvalue_max=randomwalk_params.flow_pvalue_max+1E-6, # add epsilon to workaround rounding errors
            sinks=bait_sinks,
            flow_metrics = @isdefined(tree_stats_df) ? filter(r -> r.bait_id == bait_id, tree_stats_df) : nothing,
            step_sinks=nothing,
            orig_diedges=FrameUtils.dropcategoricals(orig_diedges),
            keep_orig_diedges=false, keep_sourcesink_connected=true,
            extra_node_attrs=[:apms_p_value, :apms_median_log2, :apms_median_log2_signif,
                              :oeproteome_p_value, :oeproteome_median_log2, :oeproteome_median_log2_signif],
            kwargs...)
end

includet(joinpath(misc_scripts_path, "forceatlas3_layout.jl"))
FA = ForceAtlas3

flows_path = joinpath(analysis_path, "networks", "apms_oeproteome_flows_$(proj_info.hotnet_ver)_relaxed_prior")
isdir(flows_path) || mkdir(flows_path)
CSV.write(joinpath(flows_path, "apms_oeproteome_flows_$(proj_info.hotnet_ver)_cut_thresholds.txt"), bait_cut_thresholds_df, delim='\t')
sel_bait_ids = ["VZV-12"]
sel_valtypes = ["opt"]
sel_bait_ids = filter(r -> r.nflows == 0, innerjoin(unique!(select(bait_cut_thresholds_df, [:bait_id, :value_type, :threshold])),
                                                    tree_stats_df, on=[:bait_id, :threshold])).bait_id |> unique
sel_bait_ids = unique(first.(collect(keys(bait2cut_threshold))))
sel_valtypes = unique(last.(collect(keys(bait2cut_threshold))))
flowgraph_todo = vec([(bait_id, valtype) for bait_id in sel_bait_ids, valtype in sel_valtypes])
Threads.@threads for i in 1:length(flowgraph_todo)
    bait_id, valtype = flowgraph_todo[i]
    @info "Flowgraph #$i: bait=$bait_id cut=$valtype"
    cut_threshold_df = get(bait2cut_threshold, (bait_id, valtype), nothing)
    if isnothing(cut_threshold_df)
        @warn "  No threshold found, skipping"
        continue
    end
    cut_threshold = cut_threshold_df.threshold[1]
    flow_graph = flowgraphml(bait_id, cut_threshold,
        component_groups=false,
        flowpaths = :skip,# :flowattr,
        step_threshold=cut_threshold, maxsteps=2,
        layout_cut_threshold=cut_threshold * 1.25,
        pvalue_mw_max = 0.001, verbose=true)
    isempty(flow_graph.nodes) || HotnetUtils.layout_flowgraph!(flow_graph, scale=80, progressbar=false)
    valtype_flows_path = joinpath(flows_path, valtype)
    isdir(valtype_flows_path) || mkdir(valtype_flows_path)
    open(joinpath(valtype_flows_path, "$(proj_info.id)_$(bait_id)_flow_$(proj_info.hotnet_ver)_$(valtype).graphml"), "w") do io
        write(io, flow_graph)
    end
end

include(joinpath(misc_scripts_path, "plots", "plotly_utils.jl"))

metrics_plot_path = joinpath(flows_path, "metrics")
isdir(metrics_plot_path) || mkdir(metrics_plot_path)
for bait_id in sel_bait_ids, (metric, stat) in [(:flow_avginvlen, "max"), (:flow_avghopweight, "max"),
                                                (:flow_avgweight, "max"), (:flow_avgminedgeweight, "max"),
                                                (:maxcomponent_size, "max"), (:topn_components_sizesum, "max")]
    @info "Plotting $metric stats for $bait_id"
    realstats_df = filter(r -> !r.is_permuted && (r.bait_id == bait_id) && !ismissing(r[metric]), tree_stats_df)
    aggstats_df = filter(r -> r.bait_id == bait_id, tree_perm_aggstats_wide_df)
    isempty(aggstats_df) && continue
    extremes_df = filter(r -> (r.bait_id == bait_id) && (r.stat == stat) && (coalesce(r.quantile, NaN) == (stat == "max" ? 1.0 - sel_quantile : sel_quantile)) &&
                    (r.metric == metric) && !ismissing(r.value) && (r.value_type == "opt"), tree_extremes_df)
    bait_stats_plot = PlotlyUtils.permstats_plot(
        realstats_df, aggstats_df, extremes_df,
        threshold_range=:auto,#(max(0.0, cut_threshold_range[1]*0.75), cut_threshold_range[2]*1.25),
        yaxis_metric=metric, yaxis_log=false)
    isnothing(bait_stats_plot) && continue
    #bait_stats_plot.plot.layout[:width] = "300px"
    #bait_stats_plot.plot.layout[:height] = "200px"
    #bait_stats_plot.plot.layout[:font_size] = "20"
    try
        savefig(bait_stats_plot.plot,
            #joinpath(metrics_plot_path, "$(proj_info.id)_$(bait_id)_$(metric)_$(proj_info.hotnet_ver).pdf"))
            joinpath(metrics_plot_path, "$(bait_id)_$(metric)_linear.pdf"))
    catch e
        if e isa InterruptException
            rethrow(e)
        else
            @warn e
        end
    end
end

using Profile
using ProfileView
Profile.init(10^8, 1E-3)
pools = HHN.ObjectPools()
@time treestats_df = HHN.treecut_stats(
    bait2reactomefi_tree[sel_bait_id],
    nothing,
    walkmatrix=sel_walkmatrix,
    sources=findall(>(1.5), bait2weights[sel_bait_id].vertex),
    sinks=bait2apms_vertices[sel_bait_id], pools=pools)

ProfileView.view()

sum(length, reactomefi_tree_optcut_vzvorf[sel_bait_id])

filter(r -> r.is_hit, vertex_stats_df)
sum(vertex_stats_df.n_less .>= 90)

sel_bait_id = "VZV-4"
conncomp_stats_df = HHN.conncomponents_stats(HHN.cut(bait2reactomefi_tree[sel_bait_id], 0.005),
                                             filter(r -> r.vzvorf == sel_bait_id, vertex_stats_df),
                                             average_weights=false, mannwhitney_tests=true)
conncomp_stats_df
conncomp_stats_df = HHN.conncomponents_stats(reactomefi_tree_optcut_vzvorf[sel_bait_id],
                                             filter(r -> r.vzvorf == sel_bait_id, vertex_stats_df))

vertex_stats_df.walkpermweight_median

reactomefi_tree_optcut_stats_df = DataFrame(
    component_id = eachindex(reactomefi_tree_optcut),
    nvertices = length.(reactomefi_tree_optcut),
    nhits = sum.(v -> vertex_stats_df.is_hit[v], reactomefi_tree_optcut)
)
reactomefi_tree_optcut_stats_df.p_value = pvalue.(HypothesisTests.FisherExactTest.(
    reactomefi_tree_optcut_stats_df.nhits, sum(vertex_stats_df.is_hit),
    reactomefi_tree_optcut_stats_df.nvertices, nrow(vertex_stats_df)))
filter(r -> r.p_value <= 0.05, reactomefi_tree_optcut_stats_df)
