proj_info = (id = "vgirault_vzvapms",
             hotnet_ver = "20210707")
@info "Assembling HHotNet permuted tree results for $(proj_info.id), ver=$(proj_info.hotnet_ver)"

using Pkg
Pkg.activate(joinpath(base_scripts_path, "adhoc", proj_info.id))

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(base_scratch_path, proj_info.id)
const plots_path = joinpath(analysis_path, "plots")

using JLD2, DataFrames, CategoricalArrays, CSV, Serialization, CodecZstd, Base.Filesystem
using LinearAlgebra, Graphs, HierarchicalHotNet
using Dates
HHN = HierarchicalHotNet

vertex2gene, reactomefi_genes, _, bait2nperms, randomwalk_params =
open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver).jlser.zst"), "r") do io
    deserialize(io)
end;

daysold(filename) = Dates.value(round(Dates.now() - Dates.unix2datetime(mtime(filename)), Dates.Day(1)))

for bait_ix in eachindex(bait2nperms)
expected_bait_id, nperms = bait2nperms[bait_ix]
@info "Loading permutations for bait #$(bait_ix) (bait_id=$(expected_bait_id), nperms=$(nperms))"
bait_id, bait_weights, perm_bait_weights, sink_ixs, orig_stepmtx, orig_dwstepmtx, _, diedge_ixs, orig_vertex_walkweights, _ =
open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver)",
                                      "bait_$(bait_ix)_perms.jlser.zst"), "r") do io
    deserialize(io)
end;
@assert bait_id == expected_bait_id "Bait #$(bait_ix) is $(expected_bait_id), got data for $(bait_id)"

diedge_stepweights = vec(orig_stepmtx)[diedge_ixs]
diedge_walkweights = vec(orig_dwstepmtx)[diedge_ixs]

perm_loaded = falses(nperms)
perm_vertex_weights = perm_bait_weights.vertex
perm_vertex_walkweights = Array{Float64}(undef, size(perm_vertex_weights))
perm_diedge_stepweights = Array{Float64}(undef, (length(diedge_ixs), nperms))
perm_diedge_walkweights = similar(perm_diedge_stepweights)
perm_trees = Vector{HHN.SCCTree{Float64}}(undef, nperms)
treestats_lock = Threads.SpinLock()
perm_treestats_dfs = Vector{DataFrame}()

## assemble chunks that were calculated on the cluster using hotnet_perm_chunk.jl script
chunk_prefix = "$(proj_info.id)_hotnet_perm_$(proj_info.hotnet_ver)"
for (root, dirs, files) in Filesystem.walkdir(joinpath(scratch_path, chunk_prefix))
    bait_files = filter(file -> startswith(file, chunk_prefix * "_$(bait_ix)_") && endswith(file, ".jlser.zst") &&
                                daysold(joinpath(root, file)) <= 5, files)
    if isempty(bait_files)
        @warn "Found no bait #$(bait_ix) permutation files in $root"
        continue # skip the folder
    else
        @info "Found $(length(bait_files)) permutation file(s) for bait #$(bait_ix) in $root, processing..."
    end
    Threads.@threads for file in bait_files
        @info "Processing file $file (for bait #$bait_ix)..."
        chunk_job_info, chunk_bait_id, chunk_treeixs, chunk_vertex_walkweights,
        chunk_diedge_stepweights, chunk_diedge_walkweights,
        chunk_trees, chunk_treestats_df =
        open(ZstdDecompressorStream, joinpath(root, file), "r") do io
            deserialize(io)
        end
        @assert chunk_bait_id == bait_id "Chunk $(file) is for $(chunk_bait_id), not $(bait_id)"
        perm_loaded[chunk_treeixs] .= true
        perm_vertex_walkweights[:, chunk_treeixs] .= chunk_vertex_walkweights
        perm_diedge_stepweights[:, chunk_treeixs] .= chunk_diedge_stepweights
        perm_diedge_walkweights[:, chunk_treeixs] .= chunk_diedge_walkweights
        perm_trees[chunk_treeixs] .= chunk_trees
        lock(treestats_lock)
        push!(perm_treestats_dfs, chunk_treestats_df)
        unlock(treestats_lock)
    end
end
nperms_loaded = count(perm_loaded)
if nperms_loaded > 0
    if nperms_loaded < nperms
        @warn "$nperms_loaded of $nperms permuted trees loaded for bait #$(bait_ix) (bait_id=$(expected_bait_id))"
    else
        #@info "All of $(length(perm_loaded)) permuted trees loaded"
    end
else
    @warn "No perm chunk files for bait #$bait_ix ($bait_id) found"
    continue
end

perm_vertex_weights = perm_vertex_weights[:, perm_loaded]
perm_vertex_walkweights = perm_vertex_walkweights[:, perm_loaded]
perm_diedge_stepweights = perm_diedge_stepweights[:, perm_loaded]
perm_diedge_walkweights = perm_diedge_walkweights[:, perm_loaded]

GC.gc()

@info "Calculating diedge statistics for bait #$bait_ix ($bait_id), nde=$(length(diedge_ixs))"
diedge_stats_df = HHN.diedge_stats(length(bait_weights.vertex), diedge_ixs,
                                diedge_stepweights, diedge_walkweights,
                                perm_diedge_stepweights, perm_diedge_walkweights)
diedge_stats_df.prob_perm_walkweight_greater = 1.0 .- diedge_stats_df.ngreater_walkweight ./ nperms_loaded;
diedge_stats_df[!, :bait_id] .= bait_id;
diedge_stats_df.src_gene_id = vertex2gene[diedge_stats_df.src]
diedge_stats_df.dest_gene_id = vertex2gene[diedge_stats_df.dest]
diedge_stats_df.src_gene_name = levels!(categorical(getindex.(Ref(reactomefi_genes),
                                                              diedge_stats_df.src_gene_id)),
                                        levels(reactomefi_genes))
diedge_stats_df.dest_gene_name = levels!(categorical(getindex.(Ref(reactomefi_genes),
                                            diedge_stats_df.dest_gene_id)),
                                         levels(reactomefi_genes))
# purge diedges weights ASAP
perm_diedge_stepweights = nothing
perm_diedge_walkweights = nothing
GC.gc()

@info "Assembling tree statistics for bait #$bait_ix ($bait_id)"
perm_treestats_df = vcat(perm_treestats_dfs...)
perm_treestats_df[!, :bait_id] .= bait_id;

@info "Calculating vertex statistics for bait #$bait_ix ($bait_id), nv=$(length(bait_weights.vertex))"
vertex_stats_df = HHN.vertex_stats(bait_weights.vertex, orig_vertex_walkweights,
                            perm_vertex_weights, perm_vertex_walkweights)
vertex_stats_df.prob_perm_walkweight_greater = 1.0 .- vertex_stats_df.ngreater_walkweight ./ nperms_loaded;
vertex_stats_df[!, :bait_id] .= bait_id;
vertex_stats_df.gene_id = vertex2gene[vertex_stats_df.vertex];
vertex_stats_df.gene_name = levels!(categorical(getindex.(Ref(reactomefi_genes),
                                                           vertex_stats_df.gene_id)),
                                     levels(reactomefi_genes));

output_prefix = "$(proj_info.id)_hotnet_perm_assembled_$(proj_info.hotnet_ver)"
isdir(joinpath(scratch_path, output_prefix)) || mkdir(joinpath(scratch_path, output_prefix))

open(ZstdCompressorStream, joinpath(scratch_path, output_prefix, "$(output_prefix)_$(bait_ix).jlser.zst"), "w") do io
    serialize(io, (perm_treestats_df, vertex_stats_df, diedge_stats_df))
end

end
