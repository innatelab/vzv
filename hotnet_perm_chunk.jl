#=
job_info = (project = "vgirault_vzvapms",
            name = "vzv_hotnet_perm",
            hotnet_ver = "20210707",
            id = "vzv_hotnet_perm_20210707",
            ntrees_perchunk = 20,
            chunk = 23)
using Pkg
Pkg.activate(joinpath(base_scripts_path, "adhoc", job_info.project))
=#
job_info = (project = ARGS[1],
            name = ARGS[2],
            hotnet_ver = ARGS[3],
            id = ARGS[4],
            ntrees_perchunk = parse(Int, ARGS[5]),
            chunk = parse(Int, ARGS[6]))
@info "Permuted tree for $(job_info.project) (job '$(job_info.name)' id=$(job_info.id))" *
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
using LinearAlgebra, SparseArrays, Graphs, HierarchicalHotNet
HHN = HierarchicalHotNet

_, _, reactomefi_digraph_rev, bait2nperms, randomwalk_params =
open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver).jlser.zst"), "r") do io
    deserialize(io)
end;

bait_ids = first.(bait2nperms)
nbaittrees = cumsum(last.(bait2nperms))
chunk_1st_tree = (job_info.chunk-1) * job_info.ntrees_perchunk + 1
bait_ix = searchsortedfirst(nbaittrees, chunk_1st_tree)
@assert (bait_ix <= length(bait_ids)) "Chunk #$(job_info.chunk) outside of permuted trees range"

chunk_prefix = "$(proj_info.id)_hotnet_perm_$(proj_info.hotnet_ver)"
if isfile(joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(bait_ix)_$(job_info.chunk).jlser.zst"))
    @info "Chunk $(chunk_prefix)-#$(job_info.chunk) already exists, skipping"
    exit()
end

bait_prev_lastree = bait_ix > 1 ? nbaittrees[bait_ix-1] : 0
# indices of chunk trees to process (within bait_id perm_trees vector)
chunk_treeixs = (chunk_1st_tree - bait_prev_lastree):min(
            chunk_1st_tree - bait_prev_lastree - 1 + job_info.ntrees_perchunk,
            bait2nperms[bait_ix][2])
@info "Processing bait #$(bait_ix) ($(bait_ids[bait_ix])) permuted trees $(chunk_treeixs)"

reactomefi_mtx = LightGraphs.adjacency_matrix(reactomefi_digraph_rev, dir=:in);

bait_id, _, perm_weights, sink_ixs, _, _, _, diedge_ixs, _, _ =
open(ZstdDecompressorStream, joinpath(scratch_path, "$(proj_info.id)_hotnet_perm_input_$(proj_info.hotnet_ver)",
                                      "bait_$(bait_ix)_perms.jlser.zst"), "r") do io
    deserialize(io)
end;
@assert bait_id == bait_ids[bait_ix] "Bait #$(bait_ix) is $(bait_ids[bait_ix]), got data for $(bait_id)"

vertex_walkweights = Matrix{Float64}(undef, (size(perm_weights.vertex, 1), length(chunk_treeixs)))
diedge_stepweights = Matrix{Float64}(undef, (length(diedge_ixs), length(chunk_treeixs)))
diedge_walkweights = Matrix{Float64}(undef, (length(diedge_ixs), length(chunk_treeixs)))
trees = Vector{HHN.SCCTree{Float64}}(undef, length(chunk_treeixs))
treestats_dfs = Vector{DataFrame}()
W = HHN.weighttype(eltype(trees))
#seedling = HHN.SCCSeedling(reactomefi_mtx)

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

for (i, permix) in enumerate(chunk_treeixs)
    @info "Growing permuted tree $bait_id-#$permix ($i of $(length(chunk_treeixs)))..."
    vertex_weights = convert(Vector{W}, view(perm_weights.vertex, :, permix))
    stepmtx = HHN.stepmatrix(reactomefi_mtx, inedge_weights=view(perm_weights.edge, :, permix))
    walkmtx = HHN.random_walk_matrix(stepmtx, randomwalk_params.restart_prob)

    diff_vertex_weights = walkmtx * vertex_weights
    k = sum(diff_vertex_weights) > 0 ? sum(vertex_weights)/sum(diff_vertex_weights) : 0.0
    diff_vertex_weights .*= k
    dwstepmtx = HHN.stabilized_stepmatrix(sparse(stepmtx), vertex_weights, randomwalk_params.restart_prob,
                                          diff_vertex_weights)
    dwstepmtx = sparseclamp(dwstepmtx, absminval=randomwalk_params.step_weight_min)
    diedge_stepweights[:, i] .= vec(stepmtx)[diedge_ixs]
    diedge_walkweights[:, i] .= vec(dwstepmtx)[diedge_ixs]
    vertex_walkweights[:, i] .= diff_vertex_weights
    trees[i] = HHN.scctree(dwstepmtx, #seedling=seedling,
                           verbose=false, method=:bisect)
    @info "Treecut statistics for $bait_id-#$permix"
    spwalkmtx = sparseclamp(walkmtx, absminval=randomwalk_params.walk_weight_min)
    local treestats_df = HHN.treecut_stats(trees[i], walkmatrix=dwstepmtx, maxweight=randomwalk_params.flow_weight_max,
                                    sources=findall(>=(randomwalk_params.source_weight_min), vertex_weights),
                                    sinks=sink_ixs, sourcesinkweights=spwalkmtx, pools=nothing)
    treestats_df[!, :tree] .= permix
    treestats_df[!, :is_permuted] .= true
    push!(treestats_dfs, treestats_df)
end
treestats_df = vcat(treestats_dfs...)

chunk_prefix = "$(proj_info.id)_hotnet_perm_$(proj_info.hotnet_ver)"
isdir(joinpath(scratch_path, chunk_prefix)) || mkdir(joinpath(scratch_path, chunk_prefix))
open(ZstdCompressorStream, joinpath(scratch_path, chunk_prefix, "$(chunk_prefix)_$(bait_ix)_$(job_info.chunk).jlser.zst"), "w") do io
    serialize(io, (job_info, bait_id, chunk_treeixs, vertex_walkweights,
                   diedge_stepweights, diedge_walkweights, trees, treestats_df))
end
@info "Permuted trees chunk #$(job_info.chunk) complete"
