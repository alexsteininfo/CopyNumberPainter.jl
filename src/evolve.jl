"""
    InitialState

The copy-number state of the tree's root.

Because the upstream leaf sampler keeps the founder as the root of a sampled tree, the
root is generally *not* the most recent common ancestor of the sample. Truncal
alterations are therefore expressed as the root's initial state rather than as an
MRCA-specific special case.

Concrete states: [`Diploid`](@ref), [`Given`](@ref), [`TruncalCNAs`](@ref).
"""
abstract type InitialState end

"""
    Diploid()

Start from a normal karyotype, with the sex mode taken from the assembly. The default.
"""
struct Diploid <: InitialState end

"""
    Given(profile)

Start from `profile`, which must be on the same assembly and sex as the simulation.
Covers the real-data case of starting from a called ancestral or consensus profile.
The profile is copied, never mutated.
"""
struct Given <: InitialState
    profile::CNProfile
    # Validated here so a malformed starting state fails at the call site that made
    # it, not somewhere inside the traversal.
    function Given(profile::CNProfile)
        check_invariants(profile)
        new(profile)
    end
end

"""
    TruncalCNAs(n; wgd = 0, mode = nothing)

Apply `n` alterations, and `wgd` whole-genome doublings, to a diploid genome before
traversal begins — the truncal (clonal) copy-number state, the copy-number analogue
of clonal point mutations. The alterations are drawn from the same model as the rest
of the tree; each doubling falls at a uniformly random position among them, so a
truncal doubling can be preceded as well as followed by truncal events. Everything is
logged against the root, so it appears in the event record like any other event.

`mode` is the doubling arithmetic (`:multiply` or `:increment`, see
[`WholeGenomeDoubling`](@ref)); `nothing` takes it from the model's WGD policy, so
truncal and later doublings agree unless you say otherwise.

This is how a doubling on the trunk is expressed: the root has no incoming edge, so
[`ScheduledWGD`](@ref) cannot name it.
"""
struct TruncalCNAs <: InitialState
    n::Int
    wgd::Int
    mode::Union{Symbol,Nothing}
    function TruncalCNAs(n::Integer; wgd::Integer = 0, mode::Union{Symbol,Nothing} = nothing)
        n >= 0 || throw(ArgumentError("TruncalCNAs needs n ≥ 0, got $n"))
        wgd >= 0 || throw(ArgumentError("TruncalCNAs needs wgd ≥ 0, got $wgd"))
        mode === nothing || _check_wgd_mode(mode)
        new(Int(n), Int(wgd), mode)
    end
end

"""
    CNAModel(; rate, target, extent, kind, wgd, viability, initial)

The complete alteration model: how many alterations per edge, where each lands, how
far it runs, whether it is a gain or a loss, where whole-genome doublings fall, which
proposals are allowed, and what the root looks like.

Each component is injected and independently replaceable, and each of the three draws
also accepts a plain function. That is deliberate: downstream inference has to *fit*
these parameters, so every one of them must be addressable and cheap to vary.

# Keyword defaults
- `rate = PerDivision(1.0)` — see [`CNARate`](@ref).
- `target = CNWeighted()` — copy-number material, length × mean copy number; see
  [`TargetDraw`](@ref).
- `extent = ExtentMixture()` — focal only; set `p_arm`/`p_chromosome` to enable
  large-scale events. See [`ExtentMixture`](@ref).
- `kind = GainLoss(0.5)` — see [`KindDraw`](@ref).
- `wgd = NoWGD()` — see [`WGDPolicy`](@ref).
- `viability = RejectAndRedraw()` — see [`ViabilityRule`](@ref).
- `initial = Diploid()` — see [`InitialState`](@ref).

`rate`, `wgd`, `viability` and `initial` are type-checked at construction. `target`,
`extent` and `kind` are not, because they also accept plain functions and callable
structs.

# Examples
```julia
model = CNAModel(
    rate = PerDivision(0.5),
    target = CNWeighted(1.0),
    extent = ExtentMixture(p_chromosome = 0.05, p_arm = 0.15),
    kind = GainLoss(0.6),
    wgd = ScheduledWGD(mrca(tree, metastatic_leaves) => 1),
    initial = TruncalCNAs(4),
)
```
"""
struct CNAModel{R,T,E,K,W,V,I}
    rate::R
    target::T
    extent::E
    kind::K
    wgd::W
    viability::V
    initial::I
end

function CNAModel(; rate = PerDivision(1.0),
                    target = CNWeighted(),
                    extent = ExtentMixture(),
                    kind = GainLoss(0.5),
                    wgd = NoWGD(),
                    viability = RejectAndRedraw(),
                    initial = Diploid())
    rate isa CNARate || throw(ArgumentError(
        "rate must be a CNARate such as PerDivision(0.5), got $(repr(rate))"))
    wgd isa WGDPolicy || throw(ArgumentError(
        "wgd must be a WGDPolicy such as NoWGD() or ScheduledWGD(i => 1), got $(repr(wgd))"))
    viability isa ViabilityRule || throw(ArgumentError(
        "viability must be a ViabilityRule such as RejectAndRedraw(), got $(repr(viability))"))
    initial isa InitialState || throw(ArgumentError(
        "initial must be an InitialState such as Diploid() or TruncalCNAs(4), got $(repr(initial))"))
    CNAModel(rate, target, extent, kind, wgd, viability, initial)
end

Base.show(io::IO, m::CNAModel) = print(io, "CNAModel(rate=", m.rate, ", wgd=",
    nameof(typeof(m.wgd)), ", viability=", nameof(typeof(m.viability)),
    ", initial=", nameof(typeof(m.initial)), ")")

"""
    LoggedEvent(node, order, event)

One alteration, recorded against the edge it fell on.

`node` identifies the edge by its child, `order` is the event's position within that
edge's sequence (1-based), and `event` is the [`CNAEvent`](@ref) itself. Root-state
alterations from [`TruncalCNAs`](@ref) are logged against the root.
"""
struct LoggedEvent
    node::Int
    order::Int
    event::CNAEvent
end

"""
    CNAEvolution

The result of [`simulate_cnas`](@ref): profiles, the event log, and the rejection
tally.

# Fields
- `tree::PhyloTree` — the input tree.
- `assembly::GenomeAssembly` — the genome the simulation ran on.
- `model` — the [`CNAModel`](@ref) used.
- `profiles::Vector{Union{CNProfile,Nothing}}` — indexed by node id. Every node by
  default; leaves only when `retain_internal = false`.
- `events::Vector{LoggedEvent}` — **complete**, always, in preorder of node and then
  by `order` within a node.
- `event_ranges::Vector{UnitRange{Int}}` — indexed by node id, the contiguous run of
  `events` logged against that node (empty for a node without events).
- `root_base::CNProfile` — the root's state *before* any root-logged event (the
  [`TruncalCNAs`](@ref) alterations).
- `rejections::Dict{Symbol,Int}` — redrawn proposals, keyed by reason: a viability
  constraint, or `:no_effect` for proposals on absent DNA. Viability keys mean the realised alteration distribution is conditioned on
  viability; see [`ViabilityRule`](@ref).
- `seed::Union{Int,Nothing}` — the seed, given or drawn; `nothing` only when an explicit `rng` was passed.
- `rng_mode::Symbol` — `:global` or `:per_node`, so a run can be repeated from its record.
- `retain_internal::Bool`.

The event log is the primitive and the profiles are a cache: `root_base` plus the log
reproduce everything, so [`replay`](@ref) needs only this result, not the model, and
nothing is lost by running with `retain_internal = false`.

Its fields are part of the public API and are covered by semantic versioning.
"""
struct CNAEvolution{M}
    tree::PhyloTree
    assembly::GenomeAssembly
    model::M
    profiles::Vector{Union{CNProfile,Nothing}}
    events::Vector{LoggedEvent}
    event_ranges::Vector{UnitRange{Int}}
    root_base::CNProfile
    rejections::Dict{Symbol,Int}
    seed::Union{Int,Nothing}
    rng_mode::Symbol
    retain_internal::Bool
end

Base.show(io::IO, r::CNAEvolution) = print(io, "CNAEvolution(",
    nnodes(r.tree), " nodes, ", length(leaves(r.tree)), " leaves, ",
    length(r.events), " events, ", rejection_count(r), " rejections)")

# SplitMix64 (Steele, Lea & Flood 2014): a fixed integer mixer, so a stream derived
# from it is the same on every Julia version. `hash` gives no such promise.
function _splitmix64(x::UInt64)
    z = x + 0x9e3779b97f4a7c15
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end

# A Xoshiro stream whose four state words come straight from SplitMix64. Xoshiro(seed)
# routes the seed through Random's own seeding hash, which has changed between Julia
# versions; setting the state directly leaves only the Xoshiro256++ algorithm itself.
function _stable_rng(words::UInt64...)
    x = foldl((acc, w) -> _splitmix64(acc ⊻ w), words; init = 0x5eedc0de5eedc0de)
    g = 0x9e3779b97f4a7c15
    return Random.Xoshiro(_splitmix64(x), _splitmix64(x + g),
                          _splitmix64(x + 2g), _splitmix64(x + 3g))
end

_u64(x::Integer) = reinterpret(UInt64, Int64(x))

"""
    simulate_cnas(tree, assembly, model; rng, seed, retain_internal = true, rng_mode = :global)
        -> CNAEvolution

Draw copy-number alterations along `tree` and return every node's profile plus the
complete event log.

Traversal is depth-first and iterative, so tree depth costs heap rather than stack.
Memory is dominated by stored profiles: every node's by default, the leaves' with
`retain_internal = false`. In a lean run, the working copies along the current
root-to-node path are the only other profiles alive. Per edge, in
order: the number of doublings (from the policy's schedule, or from the edge's own stream for a `RateWGD` under `rng_mode = :per_node`), then `n_cnas(model.rate, …)`,
then the position of each doubling among the edge's events (uniformly at random), then
the segmental alterations in slot order, each drawn as target → extent → kind,
checked against the viability rule and redrawn on rejection. The realised order is
recorded, so "gained then doubled" is always distinguishable from "doubled then
gained".

Alterations here are **neutral by construction**. The tree is an input that already
encodes whatever selection produced it, so this function never kills or reweights a
cell.

# Arguments
- `tree::PhyloTree` — from [`read_newick`](@ref) or converted from a
  `NonMarkovEvolution.jl` lineage tree.
- `assembly::GenomeAssembly` — e.g. `hg38(:female)`.
- `model::CNAModel`.

# Keywords
- `rng` — an explicit random number generator. The run is then reproducible only by
  re-supplying an equal generator, and `seed` is recorded as `nothing`. Ignored if
  `seed` is given.
- `seed::Union{Int,Nothing}` — recorded in the result. When neither `seed` nor `rng` is
  given, a seed is drawn from `Random.default_rng()`, recorded and used.
- `retain_internal::Bool = true` — keep internal-node profiles. `false` keeps only
  leaves, for very large trees; the event log stays complete either way.
- `rng_mode::Symbol = :global` — `:global` threads one stream through the whole
  traversal. `:per_node` gives each edge its own stream derived from `seed` and the
  node's `source_id` (falling back to its dense id), so an edge draws identically no
  matter which other edges exist. That makes simulating on a sampled tree give
  *exactly* the same alterations as simulating on the full tree and subsetting, rather
  than only the same distribution. Requires `seed` (one is drawn if neither `seed` nor `rng` is given; an explicit `rng` is an error). A `RateWGD` then draws each
  edge's doublings from that edge's own stream; `ExactlyNWGD` cannot commute and
  warns; `ScheduledWGD(...; by = :source_id, allow_missing = true)` does. The full
  list of conditions is in the output documentation.

Streams are version-stable. A seed sets Xoshiro's state through a fixed integer mixer,
not through Julia's seeding hash, so the same seed gives the same random numbers on
every Julia version. What can still change between versions is how a *sampler* turns
those numbers into draws (for example, Distributions' `Poisson`), so a bit-identical
rerun also needs the same package versions. `save_simulation` records them.

# Examples
```julia
tree  = read_newick("lineage.nwk"; branchlength = :divisions)
res   = simulate_cnas(tree, hg38(:female), CNAModel(rate = PerDivision(0.4)); seed = 1)
grid  = BinGrid(hg38(:female), 500_000)
mat   = CNMatrix(res, grid)
write_medicc2("cells.tsv", mat)
```
"""
function simulate_cnas(tree::PhyloTree, assembly::GenomeAssembly, model::CNAModel;
                       rng::Union{Random.AbstractRNG,Nothing} = nothing,
                       seed::Union{Integer,Nothing} = nothing,
                       retain_internal::Bool = true,
                       rng_mode::Symbol = :global)
    rng_mode in (:global, :per_node) ||
        throw(ArgumentError("rng_mode must be :global or :per_node, got :$rng_mode"))
    if seed === nothing && rng === nothing
        # Drawn, recorded and used, so every run can be repeated from its own result.
        # Drawn from the default rng, so Random.seed! still controls it.
        seed = rand(Random.default_rng(), 0:typemax(Int))
    end
    seed !== nothing && !(typemin(Int) <= seed <= typemax(Int)) && throw(ArgumentError(
        "seed must fit in an Int (at most $(typemax(Int))), got $seed"))
    if seed !== nothing
        rng = _stable_rng(_u64(seed))
    elseif rng_mode === :per_node
        throw(ArgumentError(
            "rng_mode = :per_node derives each edge's stream from the seed, so pass seed " *
            "rather than an explicit rng"))
    end
    iseed = seed === nothing ? nothing : Int(seed)

    rejections = Dict{Symbol,Int}()
    events = LoggedEvent[]
    profiles = Vector{Union{CNProfile,Nothing}}(nothing, nnodes(tree))
    # Under :per_node a RateWGD draws each edge's doublings from that edge's own
    # stream, so doublings commute with sampling like everything else. ExactlyNWGD
    # cannot: which n edges it picks depends on the whole tree at hand.
    per_edge_wgd = rng_mode === :per_node && model.wgd isa RateWGD
    rng_mode === :per_node && model.wgd isa ExactlyNWGD &&
        @warn "ExactlyNWGD chooses its edges among those of this particular tree, so under " *
              "rng_mode = :per_node the doublings do not commute with leaf sampling. Use " *
              "RateWGD, or ScheduledWGD(...; by = :source_id), if they must."
    schedule = per_edge_wgd ? nothing : prepare_wgd(model.wgd, tree, rng)

    r = treeroot(tree)
    rootrng = _edge_rng(rng_mode, rng, iseed, tree, r)
    root_base = initial_profile(model.initial, assembly)
    rootp = _root_events!(copy(root_base), model.initial, model, rootrng, rejections, events, r)
    (retain_internal || isleaf(tree, r)) && (profiles[r] = rootp)

    # Iterative depth-first descent. A frame holds a node, its profile, and the index
    # of the next child to visit, so only the current root-to-node path is in memory
    # and a 10^5-deep tree cannot overflow the stack.
    fnode = Int[r]
    fprof = CNProfile[rootp]
    fnext = Int[1]
    while !isempty(fnode)
        i = fnode[end]
        kids = childrenof(tree, i)
        k = fnext[end]
        if k > length(kids)
            pop!(fnode); pop!(fprof); pop!(fnext)
            continue
        end
        fnext[end] = k + 1
        child = kids[k]
        # The parent's profile is needed again only for its later children, and is
        # stored only when retained; so the last child of an unstored parent takes it
        # over instead of copying it.
        last_use = k == length(kids) && !retain_internal
        cp = last_use ? fprof[end] : copy(fprof[end])
        erng = _edge_rng(rng_mode, rng, iseed, tree, child)
        _evolve_edge!(cp, tree, child, model, schedule, erng, rejections, events)
        (retain_internal || isleaf(tree, child)) && (profiles[child] = cp)
        push!(fnode, child); push!(fprof, cp); push!(fnext, 1)
    end

    return CNAEvolution(tree, assembly, model, profiles, events,
                        _event_ranges(events, nnodes(tree)), root_base, rejections,
                        iseed, rng_mode, retain_internal)
end

# Events are logged in preorder, grouped by node, so each node's events are one
# contiguous run of the log; this finds every run in a single pass. The log can also
# come from a file, so the layout replay and events_on rely on is checked rather than
# assumed: nodes within 1:n, one run per node, and `order` counting 1, 2, ... in a run.
function _event_ranges(events::Vector{LoggedEvent}, n::Int)
    ranges = fill(1:0, n)
    k = 1
    while k <= length(events)
        nd = events[k].node
        1 <= nd <= n || throw(ArgumentError(
            "event log refers to node $nd, but the tree has nodes 1:$n"))
        isempty(ranges[nd]) || throw(ArgumentError(
            "event log is not grouped by node: the events of node $nd are not contiguous"))
        j = k
        while j < length(events) && events[j + 1].node == nd
            j += 1
        end
        for (pos, idx) in enumerate(k:j)
            events[idx].order == pos || throw(ArgumentError(
                "events of node $nd must have order 1, 2, ... in log order; " *
                "found order $(events[idx].order) at position $pos"))
        end
        ranges[nd] = k:j
        k = j + 1
    end
    return ranges
end

function _edge_rng(mode::Symbol, rng::Random.AbstractRNG, seed::Union{Int,Nothing},
                   t::PhyloTree, i::Integer)
    mode === :global && return rng
    sid = node(t, i).source_id
    # Tagged so a node without a source_id (keyed by its dense id) cannot share a
    # stream with a node whose source_id happens to equal that dense id.
    return sid === nothing ? _stable_rng(_u64(seed), UInt64(2), _u64(i)) :
                             _stable_rng(_u64(seed), UInt64(1), _u64(sid))
end

"""
    initial_profile(init, assembly) -> CNProfile

The root's copy-number state **before** any alteration logged against the root.

This is the whole interface of an [`InitialState`](@ref): define a new subtype and a
method of this function, and both [`simulate_cnas`](@ref) and [`replay`](@ref) use it.
For [`TruncalCNAs`](@ref) it is the diploid genome the truncal events are applied to.
"""
initial_profile(::Diploid, a::GenomeAssembly) = diploid(a)
initial_profile(::TruncalCNAs, a::GenomeAssembly) = diploid(a)

function initial_profile(init::Given, a::GenomeAssembly)
    same_assembly(init.profile.assembly, a) || throw(ArgumentError(
        "Given initial profile is on $(init.profile.assembly.name)/:$(init.profile.assembly.sex) " *
        "but the simulation runs on $(a.name)/:$(a.sex)"))
    return copy(init.profile)
end

# Apply and log the alterations an initial state places at the root. Only
# TruncalCNAs has any; every other state is fully described by initial_profile.
_root_events!(p::CNProfile, ::InitialState, model, rng, rejections, events, r) = p

function _root_events!(p::CNProfile, init::TruncalCNAs, model, rng, rejections, events, r)
    wgd = WholeGenomeDoubling(something(init.mode, wgd_mode(model.wgd)))
    for (order, isdoubling) in enumerate(_wgd_positions(init.wgd, init.n, rng))
        ev = isdoubling ? wgd : _draw_cna(model, p, rng, rejections)
        apply!(p, ev)
        push!(events, LoggedEvent(r, order, ev))
    end
    return p
end

function _evolve_edge!(p::CNProfile, t::PhyloTree, i::Int, model::CNAModel,
                       schedule::Union{Dict{Int,Int},Nothing}, rng::Random.AbstractRNG,
                       rejections::Dict{Symbol,Int}, events::Vector{LoggedEvent})
    ndoublings = schedule === nothing ? n_cnas(model.wgd.rate, t, i, rng) :
                 get(schedule, i, 0)     # a read, never an iteration
    ncnas = n_cnas(model.rate, t, i, rng)
    wgd = ndoublings > 0 ? WholeGenomeDoubling(wgd_mode(model.wgd)) : nothing
    for (order, isdoubling) in enumerate(_wgd_positions(ndoublings, ncnas, rng))
        ev = isdoubling ? wgd : _draw_cna(model, p, rng, rejections)
        apply!(p, ev)
        push!(events, LoggedEvent(i, order, ev))
    end
    return p
end

# Which of an edge's `w + k` event slots hold its `w` doublings, uniformly at random.
# Placing doublings first would make "gained, then doubled" impossible within an edge,
# which matters on long edges — a trunk, or a newick edge of many divisions — where
# under :multiply a pre-doubling gain ends up +2 and a post-doubling one +1.
# The rng is consumed only when both kinds occur, so an edge without a doubling draws
# exactly the same stream whatever the WGD policy.
function _wgd_positions(w::Int, k::Int, rng::Random.AbstractRNG)
    slots = falses(w + k)
    if w > 0 && k > 0
        slots[StatsBase.sample(rng, 1:(w + k), w; replace = false)] .= true
    else
        slots .= w > 0
    end
    return slots
end

# Whether `e` would change `p` at all. Copy number 0 is absorbing, so an event that
# touches only absent DNA — gain or loss — changes nothing, and logging it would
# count an event that never physically happened.
function _has_effect(p::CNProfile, e::SegmentalCNA)
    segs = p.segments[slot(p.assembly, e.chrom, e.haplotype)]
    i = segment_index(segs, e.start)
    while i <= length(segs) && segs[i].start <= e.stop
        segs[i].cn > 0 && return true
        i += 1
    end
    return false
end

const _MAX_NO_EFFECT_REDRAWS = 10_000

function _draw_cna(model::CNAModel, p::CNProfile, rng::Random.AbstractRNG,
                   rejections::Dict{Symbol,Int})
    attempts = max_attempts(model.viability)
    tried = 0
    noeffect = 0
    while tried < attempts
        c, h = draw_target(model.target, p, rng)
        s, e, scale = draw_extent(model.extent, p, c, h, rng)
        δ = draw_kind(model.kind, p, c, h, s, e, rng)
        ev = SegmentalCNA(c, h, s, e, δ, scale)
        # A proposal on absent DNA is not an event: redraw it without spending a
        # viability attempt, so the logged count is the count of real alterations.
        if !_has_effect(p, ev)
            rejections[:no_effect] = get(rejections, :no_effect, 0) + 1
            (noeffect += 1) <= _MAX_NO_EFFECT_REDRAWS || error(
                "$(_MAX_NO_EFFECT_REDRAWS) consecutive proposals fell on absent DNA: the " *
                "profile has almost no copy-number material left to alter. Use a target rule " *
                "that weights by material, such as CNWeighted(), or stop losses earlier.")
            continue
        end
        tried += 1
        why = violation(model.viability, p, ev)
        why === nothing && return ev
        rejections[why] = get(rejections, why, 0) + 1
    end
    error("viability rejection exceeded max_attempts = $attempts: the model is proposing " *
          "almost nothing viable. Loosen the rule (lower min_total_cn, raise max_attempts), " *
          "raise p_gain so losses are less frequent, or reduce event sizes.")
end

"""
    profile(res, i) -> CNProfile

The copy-number profile of node `i`. Throws if it was not retained; use
[`replay`](@ref) to reconstruct all profiles from the event log.
"""
function profile(r::CNAEvolution, i::Integer)
    p = r.profiles[i]
    p === nothing && throw(ArgumentError(
        "no profile retained for node $i (the simulation ran with retain_internal = false); " *
        "call replay(res) to reconstruct every profile from the event log"))
    return p
end

"""
    leaf_profiles(res) -> Vector{CNProfile}

Profiles of the tree's leaves, in `leaves(res.tree)` order.
"""
leaf_profiles(r::CNAEvolution) = [profile(r, i) for i in leaves(r.tree)]

# Profiles of `ids`: from the cache where retained, otherwise from a single replay of
# the event log, so a lean run loses nothing a writer needs.
function _profiles_for(r::CNAEvolution, ids)
    all(i -> r.profiles[i] !== nothing, ids) && return [r.profiles[i] for i in ids]
    rp = replay(r)
    return [rp[i] for i in ids]
end

"""
    events_on(res, i) -> SubArray

Alterations that fell on the edge into node `i`, in the order they were applied, as a
view into `res.events`; `collect` it before modifying.
"""
events_on(r::CNAEvolution, i::Integer) = view(r.events, r.event_ranges[i])

"""
    events_below(res, i) -> Vector{LoggedEvent}

Alterations on every edge strictly below node `i` — the edges into `i`'s descendants —
in log order.
"""
function events_below(r::CNAEvolution, i::Integer)
    out = LoggedEvent[]
    stack = collect(childrenof(r.tree, i))
    below = Int[]
    while !isempty(stack)
        j = pop!(stack)
        push!(below, j)
        append!(stack, childrenof(r.tree, j))
    end
    for j in sort!(below; by = j -> first(r.event_ranges[j]))   # log order
        append!(out, view(r.events, r.event_ranges[j]))
    end
    return out
end

"""
    nevents(res) -> Int

Total number of logged alterations.
"""
nevents(r::CNAEvolution) = length(r.events)

"""
    rejection_count(res) -> Int

Total number of redrawn proposals, across all reasons: viability constraints, or
`:no_effect` for proposals that fell on absent DNA.
"""
rejection_count(r::CNAEvolution) = isempty(r.rejections) ? 0 : sum(values(r.rejections))

"""
    replay(res) -> Vector{CNProfile}

Reconstruct every node's profile from the root state plus the event log.

The event log is the primitive and the retained profiles are a cache, so this returns
exactly `res.profiles` wherever those were retained — which is what the replay test
asserts — and fills in the rest. Use it after a `retain_internal = false` run, or to
verify the traversal.
"""
function replay(r::CNAEvolution)
    t = r.tree
    out = Vector{CNProfile}(undef, nnodes(t))
    root = treeroot(t)
    base = copy(r.root_base)
    for le in view(r.events, r.event_ranges[root])
        apply!(base, le.event)
    end
    out[root] = base
    for i in preorder(t)
        i == root && continue
        p = copy(out[parentof(t, i)])
        for le in view(r.events, r.event_ranges[i])
            apply!(p, le.event)
        end
        out[i] = p
    end
    return out
end
