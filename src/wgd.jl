"""
    WGDPolicy

Where whole-genome doublings fall on the tree.

Every policy resolves, before traversal begins, to a schedule mapping a node id to the
number of doublings on its **incoming** edge. That resolution is
[`prepare_wgd`](@ref); the traversal then only ever *reads* the schedule, so no
`Dict` iteration ever consumes the random number generator.

Concrete policies: [`NoWGD`](@ref), [`ScheduledWGD`](@ref), [`ExactlyNWGD`](@ref),
[`RateWGD`](@ref).

All of them carry a `mode` — `:multiply` or `:increment`, see
[`WholeGenomeDoubling`](@ref) — which defaults to `:multiply`.
"""
abstract type WGDPolicy end

"""
    NoWGD()

No whole-genome doublings. Large-scale events can be switched off entirely this way
without changing the model's shape.
"""
struct NoWGD <: WGDPolicy end

"""
    ScheduledWGD(at; mode = :multiply, by = :id, allow_missing = false)
    ScheduledWGD(pairs...; mode = :multiply, by = :id, allow_missing = false)

Place doublings on exactly the edges you name: `at` maps a node id to the number of
doublings on the edge into it.

This is the policy for statements like *"one whole-genome doubling in the dataset, on
the edge into node i"*. Name the edge however is convenient — the tree helpers all
return dense ids:

```julia
ScheduledWGD(node_by_source_id(tree, 42) => 1)      # by upstream cell id
ScheduledWGD(mrca(tree, metastatic_leaves) => 1)    # a subclonal doubling
ScheduledWGD(node_by_label(tree, "cellA") => 2)     # two successive doublings
```

`by = :source_id` keys the schedule by upstream cell id instead of dense node id. Dense ids are renumbered by leaf sampling but source ids are not, so the same schedule then names the same edges on a full tree and on any sample of it — the form to use with `rng_mode = :per_node`. With `allow_missing = true`, ids absent from the tree are skipped: an edge that sampling pruned simply carries no observable doubling, exactly as if the full tree had been simulated and then subset. Name an `mrca` by source id with `node(full, mrca(full, leaves)).source_id`.

The root cannot be scheduled: it has no incoming edge. Use
`initial = TruncalCNAs(n; wgd = 1)` for a truncal doubling, or `initial = Given(profile)`,
to set the root state instead.
"""
struct ScheduledWGD <: WGDPolicy
    at::Dict{Int,Int}
    mode::Symbol
    by::Symbol
    allow_missing::Bool

    function ScheduledWGD(at::AbstractDict; mode::Symbol = :multiply, by::Symbol = :id,
                          allow_missing::Bool = false)
        _check_wgd_mode(mode)
        by in (:id, :source_id) || throw(ArgumentError("by must be :id or :source_id, got :$by"))
        allow_missing && by !== :source_id && throw(ArgumentError(
            "allow_missing applies only to by = :source_id: a dense id is never missing, only wrong"))
        d = Dict{Int,Int}()
        for (k, v) in at
            v >= 1 || throw(ArgumentError("scheduled doubling count for node $k must be ≥ 1, got $v"))
            d[Int(k)] = Int(v)
        end
        new(d, mode, by, allow_missing)
    end
end

ScheduledWGD(pairs::Pair...; mode::Symbol = :multiply, by::Symbol = :id,
             allow_missing::Bool = false) =
    ScheduledWGD(Dict(pairs...); mode = mode, by = by, allow_missing = allow_missing)

"""
    ExactlyNWGD(n; mode = :multiply)

Exactly `n` doublings, one each on `n` distinct non-root edges drawn uniformly at
random without replacement. Use this when the *number* of doublings matters but their
position should not be fixed. Throws if the tree has fewer than `n` non-root edges.
"""
struct ExactlyNWGD <: WGDPolicy
    n::Int
    mode::Symbol

    function ExactlyNWGD(n::Integer; mode::Symbol = :multiply)
        n >= 1 || throw(ArgumentError("ExactlyNWGD needs n ≥ 1, got $n; use NoWGD() for none"))
        _check_wgd_mode(mode)
        new(Int(n), mode)
    end
end

"""
    RateWGD(rate; mode = :multiply)
    RateWGD(λ::Real; mode = :multiply)

Doublings drawn per edge from a [`CNARate`](@ref) rule, so a doubling rate can be
expressed per division (`RateWGD(PerDivision(λ))`, or just `RateWGD(λ)`), per unit
real time (`RateWGD(PerTime(μ))`), or from the edge's mutation count. Reusing the rate
machinery means the per-division-versus-per-time question applies to doublings on the
same footing as to segmental events.
"""
struct RateWGD{R<:CNARate} <: WGDPolicy
    rate::R
    mode::Symbol

    function RateWGD(rate::R; mode::Symbol = :multiply) where {R<:CNARate}
        _check_wgd_mode(mode)
        new{R}(rate, mode)
    end
end

RateWGD(λ::Real; mode::Symbol = :multiply) = RateWGD(PerDivision(λ); mode = mode)

"""
    wgd_mode(policy) -> Symbol

The doubling arithmetic this policy applies, `:multiply` or `:increment`.
"""
wgd_mode(::NoWGD) = :multiply
wgd_mode(p::ScheduledWGD) = p.mode
wgd_mode(p::ExactlyNWGD) = p.mode
wgd_mode(p::RateWGD) = p.mode

"""
    prepare_wgd(policy, tree, rng) -> Dict{Int,Int}

Resolve `policy` into a schedule: node id ⇒ number of doublings on the edge into that
node. Called once, before traversal.
"""
prepare_wgd(::NoWGD, ::PhyloTree, ::Random.AbstractRNG) = Dict{Int,Int}()

function prepare_wgd(p::ScheduledWGD, t::PhyloTree, ::Random.AbstractRNG)
    out = Dict{Int,Int}()
    for k in sort!(collect(keys(p.at)))
        if p.by === :source_id
            i = get(t.source_index, k, nothing)
            if i === nothing
                p.allow_missing && continue
                throw(ArgumentError("scheduled WGD names source_id $k, which no node of this tree " *
                                    "carries; pass allow_missing = true if the edge may be absent"))
            end
        else
            1 <= k <= nnodes(t) || throw(ArgumentError(
                "scheduled WGD names node $k, which is not in this tree (1:$(nnodes(t)))"))
            i = k
        end
        isroot(t, i) && throw(ArgumentError(
            "node $i is the root and has no incoming edge; for a truncal doubling use " *
            "initial = TruncalCNAs(n; wgd = 1), or set the root state with initial = Given(profile)"))
        out[i] = p.at[k]
    end
    return out
end

function prepare_wgd(p::ExactlyNWGD, t::PhyloTree, rng::Random.AbstractRNG)
    candidates = [i for i in preorder(t) if !isroot(t, i)]
    length(candidates) >= p.n || throw(ArgumentError(
        "ExactlyNWGD($(p.n)) needs $(p.n) non-root edges but the tree has only $(length(candidates))"))
    chosen = StatsBase.sample(rng, candidates, p.n; replace = false)
    return Dict{Int,Int}(i => 1 for i in chosen)
end

function prepare_wgd(p::RateWGD, t::PhyloTree, rng::Random.AbstractRNG)
    out = Dict{Int,Int}()
    for i in preorder(t)          # deterministic order, so the schedule is reproducible
        isroot(t, i) && continue
        k = n_cnas(p.rate, t, i, rng)
        k > 0 && (out[i] = k)
    end
    return out
end
