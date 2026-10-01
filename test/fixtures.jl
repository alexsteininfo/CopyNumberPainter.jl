# Shared synthetic fixtures. No real data, ever: see the data rule in README.md.

"""
    toy_assembly(; nchrom = 2, len = 1000, sex = :female)

A tiny assembly for tests: `nchrom` chromosomes of `len` bp with a centromere in the
middle fifth. Chromosome names are `"chr1"`, `"chr2"`, … so sex handling is *not*
triggered; use `toy_sex_assembly` for that.
"""
function toy_assembly(; nchrom::Int = 2, len::Int = 1000, sex::Symbol = :female)
    specs = [CopyNumberPainter.ChromosomeSpec("chr$(i)", len,
                (2 * len ÷ 5 + 1):(3 * len ÷ 5)) for i in 1:nchrom]
    CopyNumberPainter.GenomeAssembly("toy", sex, specs, fill(2, nchrom))
end

"""
    toy_sex_assembly(sex; len = 1000)

Two autosomes plus `chrX` and `chrY`, so hemizygosity and zero-ploidy chromosomes are
exercised. `sex` is `:female` or `:male`.
"""
function toy_sex_assembly(sex::Symbol; len::Int = 1000)
    names = ["chr1", "chr2", "chrX", "chrY"]
    specs = [CopyNumberPainter.ChromosomeSpec(n, len,
                (2 * len ÷ 5 + 1):(3 * len ÷ 5)) for n in names]
    CopyNumberPainter.GenomeAssembly("toysex", sex, specs)
end

"""
    hemizygous_assembly(; len = 1000)

One chromosome present in a single copy. Used for the pathological-rejection test:
the only whole-chromosome loss available drives total copy number to zero.
"""
function hemizygous_assembly(; len::Int = 1000)
    specs = [CopyNumberPainter.ChromosomeSpec("chr1", len,
                (2 * len ÷ 5 + 1):(3 * len ÷ 5))]
    CopyNumberPainter.GenomeAssembly("hemi", :male, specs, [1])
end

"""
    induced_subtree(tree, keep_leaves) -> PhyloTree

Test-only stand-in for `NonMarkovEvolution.sample_leaves`: keep `keep_leaves` plus
every ancestor of a kept leaf, **retaining unary nodes and keeping the founder as the
root**.

This is a fixture, not package API — leaf sampling belongs upstream. It keeps the "sampling commutes" tests independent of `NonMarkovEvolution.jl`; `test_ext.jl` repeats them with the real `sample_leaves` when that package is available. The
"prune but never collapse" behaviour is the load-bearing part: every division
ancestral to a kept leaf must remain a node, or a kept cell's root-to-leaf path would
lose alteration-drawing opportunities. `source_id`s are preserved, which is what makes
edges identifiable across the two trees.
"""
function induced_subtree(t::PhyloTree, keep_leaves::AbstractVector{<:Integer})
    keepset = Set{Int}()
    for l in keep_leaves, anc in ancestors(t, l)
        push!(keepset, anc)
    end
    old = sort!(collect(keepset))
    newid = Dict(o => i for (i, o) in enumerate(old))
    parents = Vector{Union{Int,Nothing}}(undef, length(old))
    for (i, o) in enumerate(old)
        p = parentof(t, o)
        parents[i] = p === nothing ? nothing : newid[p]
    end
    return phylotree(parents;
        birthtimes     = [node(t, o).birthtime      for o in old],
        edge_divisions = [node(t, o).edge_divisions for o in old],
        edge_mutations = [node(t, o).edge_mutations for o in old],
        labels         = [node(t, o).label          for o in old],
        source_ids     = [something(node(t, o).source_id, o) for o in old])
end

"""
    binary_lineage(depth; divisions = 1, dt = 1.0) -> PhyloTree

A complete binary tree of the given `depth` (so `2^depth` leaves), with `source_id`
equal to the dense id, one division per edge, and birthtimes advancing by `dt` per
level.
"""
function binary_lineage(depth::Int; divisions::Int = 1, dt::Float64 = 1.0)
    n = 2^(depth + 1) - 1
    parents = Vector{Union{Int,Nothing}}(undef, n)
    bt = Vector{Float64}(undef, n)
    parents[1] = nothing
    bt[1] = 0.0
    for i in 2:n
        p = i ÷ 2
        parents[i] = p
        bt[i] = bt[p] + dt
    end
    return phylotree(parents;
        birthtimes = bt,
        edge_divisions = vcat(nothing, fill(divisions, n - 1)),
        edge_mutations = vcat(nothing, fill(2 * divisions, n - 1)),
        source_ids = collect(1:n))
end

"""
    ZeroRNG()

An RNG whose `rand()` is always exactly `0.0`, the one value a `u <= acc` weighted
sampler gets wrong.
"""
struct ZeroRNG <: Random.AbstractRNG end
Random.rand(::ZeroRNG, ::Random.SamplerTrivial{Random.CloseOpen01{Float64}}) = 0.0

"""
    DoubledChr1()

A user-defined `InitialState`: chromosome 1, haplotype 1 at copy number 2. It exists to
test that `initial_profile` is the whole extension interface.
"""
struct DoubledChr1 <: CopyNumberPainter.InitialState end
function CopyNumberPainter.initial_profile(::DoubledChr1, a::CopyNumberPainter.GenomeAssembly)
    p = diploid(a)
    apply!(p, SegmentalCNA(1, 1, 1, chromlength(a, 1), 1, :chromosome))
    return p
end

"""
    CyclicUserStruct()

A mutable user struct that can hold a reference to itself, used to test that
`_describe` handles cycles without stack overflow.
"""
mutable struct CyclicUserStruct
    self
end

"""
    TypeFieldStruct()

A user struct with a field holding a type (e.g. `Int`), used to test that
`_describe` handles type-valued fields without recursing into their structure.
"""
struct TypeFieldStruct
    type_field::Type
    value::Int
end

"""
    MatrixFieldStruct()

A user struct with a matrix field, used to test that `_describe` uses summary
for non-vector arrays instead of dumping all elements.
"""
struct MatrixFieldStruct
    mat::Matrix{Int}
end

"""
    LongVectorFieldStruct()

A user struct with a long vector field, used to test that `_describe` uses summary
for vectors with >1000 elements instead of dumping all elements.
"""
struct LongVectorFieldStruct
    vec::Vector{Int}
end

"""
    TupleFieldStruct()

A user struct with a Tuple field, used to test that `_describe` converts Tuples
to described vectors.
"""
struct TupleFieldStruct
    tup::Tuple
end

"""
    SetFieldStruct()

A user struct with a Set field, used to test that `_describe` converts Sets
to sorted described vectors when sortable.
"""
struct SetFieldStruct
    s::Set{Int}
end

"""
    SelfContainingVectorStruct()

A mutable user struct containing a vector that references the struct itself,
used to test that `_describe` handles cycles in container types.
"""
mutable struct SelfContainingVectorStruct
    vec::Vector{Any}
end
