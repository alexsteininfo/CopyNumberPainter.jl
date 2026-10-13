"""
    CopyNumberPainterNonMarkovEvolutionExt

Bridge from `NonMarkovEvolution.jl` (0.5) lineage trees to `CopyNumberPainter.PhyloTree`.
A `NonMarkovEvolution.CellNode` is a handle on one row of a struct-of-arrays
`LineageTree`; the converter walks the subtree below it through the public node
properties (`left`, `right`, `id`, `birthtime`, `mutations`).

This is a **package extension**, loaded only when both packages are present. Loading
`CopyNumberPainter` alone gives the copy-number modeller with no simulator in the
dependency chain — which is the point, because the downstream inference package must
be installable and runnable against real patient data with no simulator anywhere in
its dependencies.
"""
module CopyNumberPainterNonMarkovEvolutionExt

using CopyNumberPainter
using CopyNumberPainter: PhyloNode, PhyloTree, founder_mutations
using NonMarkovEvolution: CellNode, LeafSample

"""
    PhyloTree(root::CellNode) -> PhyloTree

Convert a `NonMarkovEvolution.jl` lineage tree, full or sampled: the subtree below
`root`, which becomes the root of the result. For a population use
`PhyloTree(single_root(pop))`.

Field mapping:

| `PhyloNode` field | source |
|:---|:---|
| `birthtime` | `node.birthtime` |
| `edge_divisions` | `1` — one lineage-tree edge is exactly one division |
| `edge_mutations` | `node.mutations`, the mutations acquired at this cell's birth |
| `source_id` | `node.id`, so `node_by_source_id` keeps working after leaf sampling |
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
function CopyNumberPainter.PhyloTree(root::CellNode)
    par = Union{Int,Nothing}[]
    kids = Vector{Int}[]
    bt = Union{Float64,Nothing}[]
    divs = Union{Int,Nothing}[]
    muts = Union{Int,Nothing}[]
    sids = Union{Int,Nothing}[]

    # Iterative preorder. Push the right child first so the left is popped first and
    # children end up in left-to-right order.
    stack = Tuple{CellNode,Union{Int,Nothing}}[(root, nothing)]
    while !isempty(stack)
        nd, p = pop!(stack)
        push!(par, p)
        push!(kids, Int[])
        i = length(par)
        p === nothing || push!(kids[p], i)
        # `node.data` reads all five columns at once; the children are separate lookups.
        cell = nd.data
        push!(bt, Float64(cell.birthtime))
        push!(divs, p === nothing ? nothing : 1)
        push!(muts, p === nothing ? nothing : Int(cell.mutations))
        push!(sids, Int(cell.id))
        r, l = nd.right, nd.left
        r === nothing || push!(stack, (r, i))
        l === nothing || push!(stack, (l, i))
    end

    nodes = [PhyloNode(i, par[i], kids[i], bt[i], divs[i], muts[i], nothing, sids[i])
             for i in eachindex(par)]
    return PhyloTree(nodes)
end

"""
    PhyloTree(s::LeafSample) -> PhyloTree

Convert the lineage tree of a `NonMarkovEvolution.sample_leaves` draw. Sampling
prunes without collapsing and copies the kept rows into a tree of their own, so the
result keeps one edge per real division and the founder as root; see
[`PhyloTree`](@ref).
"""
CopyNumberPainter.PhyloTree(s::LeafSample) = CopyNumberPainter.PhyloTree(s.root)

"""
    founder_mutations(root::CellNode) -> Int

Mutations the founder cell acquired at its own birth. See
`CopyNumberPainter.founder_mutations`.
"""
CopyNumberPainter.founder_mutations(root::CellNode) = Int(root.mutations)

end # module
