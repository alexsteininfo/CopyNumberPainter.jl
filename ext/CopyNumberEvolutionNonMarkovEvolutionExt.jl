"""
    CopyNumberEvolutionNonMarkovEvolutionExt

Bridge from `NonMarkovEvolution.jl`'s pointer-based lineage trees to
`CopyNumberEvolution.PhyloTree`.

This is a **package extension**, loaded only when both packages are present. Loading
`CopyNumberEvolution` alone gives the copy-number modeller with no simulator in the
dependency chain — which is the point, because the downstream inference package must
be installable and runnable against real patient data with no simulator anywhere in
its dependencies.
"""
module CopyNumberEvolutionNonMarkovEvolutionExt

using CopyNumberEvolution
using CopyNumberEvolution: PhyloNode, PhyloTree, founder_mutations
using NonMarkovEvolution: BinaryNode, NonMarkovCell, LeafSample

"""
    PhyloTree(root::BinaryNode{NonMarkovCell}) -> PhyloTree

Convert a `NonMarkovEvolution.jl` lineage tree, full or sampled.

Field mapping:

| `PhyloNode` field | source |
|:---|:---|
| `birthtime` | `cell.birthtime` |
| `edge_divisions` | `1` — one lineage-tree edge is exactly one division |
| `edge_mutations` | `cell.drivers`, the driver mutations acquired at this cell's birth |
| `source_id` | `cell.id`, so `node_by_source_id` keeps working after leaf sampling |
| `label` | `nothing` |

All three rate rules therefore apply to a converted tree.

Two things this deliberately does **not** do. It does not sample: sampling is
`NonMarkovEvolution.jl`'s own operation, and this converts whatever tree it is
handed. And it does not prune or collapse: unary nodes are preserved, because a
sampled cell's root-to-leaf path must keep one alteration-drawing opportunity per real
division. Pass a `sample_leaves` draw directly to convert a sample.

The root's `edge_mutations` is `nothing`, since the founder has no incoming edge — see
[`founder_mutations`](@ref) if you want those mutations translated into truncal
alterations.
"""
function CopyNumberEvolution.PhyloTree(root::BinaryNode{NonMarkovCell})
    par = Union{Int,Nothing}[]
    kids = Vector{Int}[]
    bt = Union{Float64,Nothing}[]
    divs = Union{Int,Nothing}[]
    muts = Union{Int,Nothing}[]
    sids = Union{Int,Nothing}[]

    # Iterative preorder. Push the right child first so the left is popped first and
    # children end up in left-to-right order.
    stack = Tuple{BinaryNode{NonMarkovCell},Union{Int,Nothing}}[(root, nothing)]
    while !isempty(stack)
        nd, p = pop!(stack)
        push!(par, p)
        push!(kids, Int[])
        i = length(par)
        p === nothing || push!(kids[p], i)
        cell = nd.data
        push!(bt, Float64(cell.birthtime))
        push!(divs, p === nothing ? nothing : 1)
        push!(muts, p === nothing ? nothing : Int(cell.drivers))
        push!(sids, Int(cell.id))
        nd.right === nothing || push!(stack, (nd.right, i))
        nd.left === nothing || push!(stack, (nd.left, i))
    end

    nodes = [PhyloNode(i, par[i], kids[i], bt[i], divs[i], muts[i], nothing, sids[i])
             for i in eachindex(par)]
    return PhyloTree(nodes)
end

"""
    PhyloTree(s::LeafSample) -> PhyloTree

Convert the lineage tree of a `NonMarkovEvolution.sample_leaves` draw. Sampling
prunes without collapsing, so the result keeps one edge per real division and the
founder as root; see [`PhyloTree`](@ref).
"""
CopyNumberEvolution.PhyloTree(s::LeafSample) = CopyNumberEvolution.PhyloTree(s.root)

"""
    founder_mutations(root::BinaryNode{NonMarkovCell}) -> Int

Driver mutations the founder cell acquired at its own birth. See
`CopyNumberEvolution.founder_mutations`.
"""
CopyNumberEvolution.founder_mutations(root::BinaryNode{NonMarkovCell}) =
    Int(root.data.drivers)

end # module
