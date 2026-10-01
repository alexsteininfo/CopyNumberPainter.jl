"""
    CopyNumberEvolutionAbstractTreesExt

The `AbstractTrees.jl` interface for `CopyNumberEvolution.NodeRef`, loaded only when
both packages are present.
"""
module CopyNumberEvolutionAbstractTreesExt

using CopyNumberEvolution
using CopyNumberEvolution: NodeRef
import AbstractTrees

AbstractTrees.children(n::NodeRef) = [NodeRef(n.tree, c) for c in childrenof(n.tree, n.id)]
AbstractTrees.nodevalue(n::NodeRef) = n.id
AbstractTrees.ParentLinks(::Type{NodeRef}) = AbstractTrees.StoredParents()
AbstractTrees.parent(n::NodeRef) =
    (p = parentof(n.tree, n.id); p === nothing ? nothing : NodeRef(n.tree, p))
AbstractTrees.NodeType(::Type{NodeRef}) = AbstractTrees.HasNodeType()
AbstractTrees.nodetype(::Type{NodeRef}) = NodeRef
AbstractTrees.printnode(io::IO, n::NodeRef) = print(io, cellname(n.tree, n.id))

end # module
