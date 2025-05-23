#=
job_info = (project = "vgirault_vzvapms",
            name = "vzv_hotnet",
            hotnet_ver = "20210707",
            id = "vzv_hotnet_20210707",
            nbaits_perchunk = 2,
            chunk = 20)
using Pkg
Pkg.activate(joinpath(base_scripts_path, "adhoc", job_info.project))
=#

job_info = (project = ARGS[1],
            name = ARGS[2],
            hotnet_ver = ARGS[3],
            id = ARGS[4],
            nbaits_perchunk = parse(Int, ARGS[5]),
            chunk = parse(Int, ARGS[6]))
@info "Treestats for $(job_info.project) (job '$(job_info.name)' id=$(job_info.id))" *
      " chunk #$(job_info.chunk)"

proj_info = (id = job_info.project,
             hotnet_ver = job_info.hotnet_ver)

const misc_scripts_path = joinpath(base_scripts_path, "misc_jl");
const analysis_path = joinpath(base_analysis_path, proj_info.id)
const data_path = joinpath(analysis_path, "data")
const results_path = joinpath(analysis_path, "results")
const scratch_path = joinpath(base_scratch_path, proj_info.id)
const plots_path = joinpath(analysis_path, "plots")

using JLD2, DataFrames, CSV, Serialization, CodecZstd
using LinearAlgebra, Graphs, HierarchicalHotNet
HHN = HierarchicalHotNet

_, _, reactomefi_digraph_rev, bait2nperms, randomwalk_params =
open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver).jlser.zst"), "r") do io
    deserialize(io)
end;

bait_ids = first.(bait2nperms)
bait_ix1 = (job_info.chunk-1) * job_info.nbaits_perchunk + 1
@assert (bait_ix1 <= length(bait_ids)) "Chunk #$(job_info.chunk): outside of valid range"
bait_ixs = bait_ix1:min(bait_ix1 + job_info.nbaits_perchunk - 1, length(bait_ids))

reactomefi_mtx = Matrix(LightGraphs.adjacency_matrix(reactomefi_digraph_rev, dir=:in));

treestats_dfs = Vector{DataFrame}()
for (i, bait_ix) in enumerate(bait_ixs)
    @info "Processing bait #$(bait_ix) ($(bait_ids[bait_ix])), $i of $(length(bait_ixs)))..."
    bait_id, bait_weights, _, sink_ixs, _, dwstepmtx, walkmtx, diedge_ixs, _, tree =
        open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver)",
                                        "bait_$(bait_ix)_perms.jlser.zst"), "r") do io
        deserialize(io)
    end
    @assert bait_id == bait_ids[bait_ix] "Bait #$(bait_ix) is $(bait_ids[bait_ix]), got data for $(bait_id)"
    @info "Treecut statistics for bait #$(bait_ix) ($(bait_id), $i of $(length(bait_ixs)))..."
    local treestats_df = HHN.treecut_stats(tree, walkmatrix=dwstepmtx, maxweight=randomwalk_params.flow_weight_max,
                                    sources=findall(>=(randomwalk_params.source_weight_min), bait_weights.vertex),
                                    sinks=sink_ixs, sourcesinkweights=walkmtx, pools=nothing)
    treestats_df[!, :tree] .= 0
    treestats_df[!, :is_permuted] .= false
    treestats_df[!, :bait_id] .= bait_id
    push!(treestats_dfs, treestats_df)
end
treestats_df = vcat(treestats_dfs...)

chunk_prefix = "$(proj_info.id)_hotnet_treestats_$(proj_info.hotnet_ver)"
isdir(joinpath(scratch_path, chunk_prefix)) || mkdir(joinpath(scratch_path, chunk_prefix))
open(ZstdCompressorStream, joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(job_info.chunk).jlser.zst"), "w") do io
    serialize(io, (job_info, bait_ids[bait_ixs], treestats_df))
end
@info "Treestats chunk #$(job_info.chunk) complete"
