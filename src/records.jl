# TOML-safe descriptions of the objects a saved simulation must record. A model can
# hold closures and user types that cannot be rebuilt from a file, so the record is
# descriptive: enough to know exactly what ran, not to run it again.

# Every method carries a depth so that a self-referencing or very deep object ends as a
# string instead of overflowing the stack; the one-argument form is the public entry.
_describe(x::Union{Bool,Integer,AbstractFloat,AbstractString}, ::Int) = x
_describe(x::Symbol, ::Int) = String(x)
_describe(::Nothing, ::Int) = "nothing"
_describe(x::Type, ::Int) = string(x)
_describe(x::UnitRange, ::Int) = [first(x), last(x)]
_describe(f::Function, ::Int) = Dict{String,Any}("type" => "function", "repr" => string(f))
_describe(a::GenomeAssembly, ::Int) = _assembly_meta(a)
# A profile is fully recorded in `<prefix>_root.tsv` and `_profiles.tsv`; dumping its
# segments here would bloat the metadata of every run that starts from a given state.
_describe(p::CNProfile, ::Int) = "CNProfile($(p.assembly.name), $(sum(length, p.segments)) segments)"

function _describe(x::AbstractVector, depth::Int)
    depth > 8 && return string(x)
    length(x) > 1000 && return summary(x)
    return Any[_describe(v, depth + 1) for v in x]
end

function _describe(x::AbstractDict, depth::Int)
    depth > 8 && return string(x)
    return Dict{String,Any}(string(k) => _describe(v, depth + 1) for (k, v) in x)
end

function _describe(x::Tuple, depth::Int)
    depth > 8 && return string(x)
    return Any[_describe(v, depth + 1) for v in x]
end

function _describe(x::AbstractSet, depth::Int)
    depth > 8 && return string(x)
    elements = collect(x)
    # Set iteration order is arbitrary, so sort when the elements allow it: the same
    # model must give the same metadata file.
    try
        sort!(elements)
    catch
        # unsortable elements keep the set's own order
    end
    return Any[_describe(v, depth + 1) for v in elements]
end

# Matrices and higher arrays are data, not configuration: record their shape only.
_describe(x::AbstractArray, ::Int) = summary(x)

function _describe(x, depth::Int)
    depth > 8 && return string(x)
    T = typeof(x)
    isstructtype(T) || return string(x)
    d = Dict{String,Any}("type" => string(nameof(T)))
    for f in fieldnames(T)
        d[string(f)] = _describe(getfield(x, f), depth + 1)
    end
    return d
end

_describe(x) = _describe(x, 0)

_assembly_meta(a::GenomeAssembly) = Dict{String,Any}(
    "name" => a.name, "sex" => String(a.sex),
    "chromosomes" => [c.name for c in a.chromosomes],
    "lengths" => [c.length for c in a.chromosomes],
    "centromere_start" => [first(c.centromere) for c in a.chromosomes],
    "centromere_stop" => [last(c.centromere) for c in a.chromosomes],
    "ploidy" => copy(a.ploidy))

_assembly_from_meta(d::AbstractDict) = GenomeAssembly(d["name"], Symbol(d["sex"]),
    [ChromosomeSpec(n, Int(l), Int(s):Int(e)) for (n, l, s, e) in
         zip(d["chromosomes"], d["lengths"], d["centromere_start"], d["centromere_stop"])],
    Int[Int(x) for x in d["ploidy"]])

"""
    ModelRecord(description, repr)

The model of a simulation read back by [`load_simulation`](@ref): its recorded
description (component types and parameters, as saved in the metadata TOML) and its
original printed form. A model can hold closures and user types that a file cannot
rebuild, so a loaded result describes the model it ran rather than containing it;
everything computed *from* the result — profiles, replay, projection, export — works
as for a fresh one.
"""
struct ModelRecord
    description::Dict{String,Any}
    repr::String
end

Base.show(io::IO, r::ModelRecord) = print(io, r.repr)
_describe(r::ModelRecord, depth::Int) = r.description

const _FORMAT_VERSION = 1

# One saved part: a missing or unreadable optional file gives a warning and `nothing`,
# so the remaining parts still load (the convention of EcDNAModelling.load_measurement).
function _tryload(read_part, path)
    isfile(path) || (@warn "saved part $path is missing; it is rebuilt or left empty"; return nothing)
    try
        return read_part(path)
    catch e
        @warn "could not read saved part $path; it is rebuilt or left empty" exception = e
        return nothing
    end
end

function _check_bundle(meta, format, prefix)
    # Bundles written before the package was renamed carry the old package name in
    # their format tag; the files themselves are unchanged, so they still load.
    legacy = replace(format, "CopyNumberPainter" => "CopyNumberEvolution")
    get(meta, "format", "") in (format, legacy) || throw(ArgumentError(
        "$(prefix)_meta.toml is not a $format bundle"))
    # Int(...): TOML hands back an Int128 for integers that do not fit in an Int64.
    Int(meta["format_version"]) <= _FORMAT_VERSION || throw(ArgumentError(
        "$(prefix)_meta.toml has format version $(meta["format_version"]), newer than the " *
        "$(_FORMAT_VERSION) this version of CopyNumberPainter reads; upgrade the package"))
    return meta
end

_version_meta() = Dict{String,Any}("format_version" => _FORMAT_VERSION,
    "package_version" => string(pkgversion(@__MODULE__)), "julia_version" => string(VERSION))

function _check_prefix(prefix)
    dir = dirname(prefix)
    (isempty(dir) || isdir(dir)) || throw(ArgumentError("directory $dir does not exist"))
end

# The first branch-length field every non-root node carries, or :none.
function _complete_branchlength(t::PhyloTree)
    nonroot = [node(t, i) for i in 1:nnodes(t) if !isroot(t, i)]
    all(nd -> nd.edge_divisions !== nothing, nonroot) && return :divisions
    all(nd -> nd.birthtime !== nothing, t.nodes) && return :time
    all(nd -> nd.edge_mutations !== nothing, nonroot) && return :mutations
    return :none
end

"""
    save_simulation(prefix, res; profiles = :all, compress = false, newick = :auto) -> prefix

Write a [`CNAEvolution`](@ref) as a set of files sharing `prefix` (the directory
must exist; existing files are overwritten):

- `<prefix>_meta.toml` — everything needed to interpret and reproduce the run: seed,
  `rng_mode`, package and Julia versions, the assembly, a description of every model
  component, the rejection tally, counts, and the list of files written (always);
- `<prefix>_tree.tsv` — the tree, losslessly ([`write_tree`](@ref)) (always);
- `<prefix>_events.tsv` — the complete event log (always);
- `<prefix>_root.tsv` — the root's state before any root-logged event (always);
- `<prefix>_profiles.tsv` — every node's profile (`profiles = :all`), the tips only
  (`:leaves`), or nothing (`:none`);
- `<prefix>.nwk` — the tree in newick, for other tools, carrying the first branch-length
  field that every edge has (`newick = :auto`), a named one, or none (`:none`).

`compress = true` gzips the event and profile tables. The event log is the
primitive, so [`load_simulation`](@ref) rebuilds any profile that was not saved.
"""
function save_simulation(prefix::AbstractString, res::CNAEvolution;
                         profiles::Symbol = :all, compress::Bool = false, newick::Symbol = :auto)
    _check_prefix(prefix)
    profiles in (:all, :leaves, :none) ||
        throw(ArgumentError("profiles must be :all, :leaves or :none, got :$profiles"))
    newick in (:auto, :none, BRANCHLENGTH_MODES...) ||
        throw(ArgumentError("newick must be :auto, :none, :divisions, :time or :mutations, got :$newick"))
    # Everything that can refuse is checked before the first file is opened, so a refusal
    # leaves no partial bundle behind.
    for i in 1:nnodes(res.tree)
        _tree_label_field(res.tree, i)
        _tsv_field(cellname(res.tree, i))
    end
    nwk = newick === :auto ? _complete_branchlength(res.tree) : newick
    nwkstring = nwk === :none ? nothing : newick_string(res.tree; branchlength = nwk)
    base = basename(prefix)
    ext = compress ? ".tsv.gz" : ".tsv"
    files = Dict{String,Any}()
    part(key, suffix) = (files[key] = base * suffix; prefix * suffix)

    write_tree(part("tree", "_tree.tsv"), res.tree)
    write_events(part("events", "_events" * ext), res)
    r = treeroot(res.tree)
    _write_profile_table(part("root", "_root.tsv"), res.tree, res.assembly, [r], [res.root_base])
    if profiles !== :none
        ids = profiles === :all ? collect(1:nnodes(res.tree)) : leaves(res.tree)
        write_profiles(part("profiles", "_profiles" * ext), res; cells = ids)
    end
    nwkstring === nothing || _with_io(io -> print(io, nwkstring), part("newick", ".nwk"))

    meta = merge(_version_meta(), Dict{String,Any}(
        "format" => "CopyNumberPainter.simulation",
        "outputs" => sort!(collect(keys(files))), "files" => files,
        "rng_mode" => String(res.rng_mode), "retain_internal" => res.retain_internal,
        "nnodes" => nnodes(res.tree), "nleaves" => length(leaves(res.tree)),
        "nevents" => nevents(res), "saved_profiles" => String(profiles),
        "newick_branchlength" => String(nwk),
        "rejections" => Dict{String,Any}(String(k) => v for (k, v) in res.rejections),
        "assembly" => _assembly_meta(res.assembly),
        "model" => _describe(res.model), "model_repr" => sprint(show, res.model)))
    res.seed === nothing || (meta["seed"] = res.seed)
    open(io -> TOML.print(io, meta; sorted = true), prefix * "_meta.toml", "w")
    return prefix
end

"""
    load_simulation(prefix) -> CNAEvolution{ModelRecord}

Read the files written by [`save_simulation`](@ref)`(prefix, res)`. The metadata, tree,
event log and root files are required. A profile table that is missing or unreadable
gives a warning; every profile not read from it is rebuilt from the event log, so
the result is always complete. The model comes back as a [`ModelRecord`](@ref).
"""
function load_simulation(prefix::AbstractString)
    meta = _check_bundle(TOML.parsefile(prefix * "_meta.toml"), "CopyNumberPainter.simulation", prefix)
    file(key) = joinpath(dirname(prefix), meta["files"][key])
    a = _assembly_from_meta(meta["assembly"])
    tree = read_tree(file("tree"))
    events = read_events(file("events"), a)
    root_base = only(values(read_profiles(file("root"), a)))
    profs = Vector{Union{CNProfile,Nothing}}(nothing, nnodes(tree))
    if haskey(meta["files"], "profiles")
        loaded = _tryload(p -> read_profiles(p, a), file("profiles"))
        if loaded !== nothing
            for (id, p) in loaded
                profs[id] = p
            end
        end
    end
    # Int(...): a drawn seed can exceed 2^62, which TOML reads back as an Int128.
    seed = haskey(meta, "seed") ? Int(meta["seed"]) : nothing
    res = CNAEvolution(tree, a, ModelRecord(meta["model"], meta["model_repr"]), profs, events,
                       _event_ranges(events, nnodes(tree)), root_base,
                       Dict{Symbol,Int}(Symbol(k) => Int(v) for (k, v) in meta["rejections"]),
                       seed, Symbol(meta["rng_mode"]), true)
    if any(isnothing, profs)                  # the log is the primitive: fill in the rest
        rp = replay(res)
        for i in eachindex(profs)
            profs[i] === nothing && (profs[i] = rp[i])
        end
    end
    return res
end

const _BIN_RULES = Dict("LengthWeightedMajority" => LengthWeightedMajority(),
                        "AreaWeightedMean" => AreaWeightedMean())

"""
    save_matrix(prefix, m; compress = false) -> prefix

Write a [`CNMatrix`](@ref) as files sharing `prefix`: `<prefix>_meta.toml` (the
assembly, grid size and mask, bin rule, phasing, versions and file list),
`<prefix>_bins.tsv` ([`write_bins`](@ref)), `<prefix>_total.tsv` and, for an
allele-specific matrix, `<prefix>_allele_<h>.tsv` per haplotype track
([`write_matrix`](@ref)). `compress = true` gzips the matrix tables.
"""
function save_matrix(prefix::AbstractString, m::CNMatrix; compress::Bool = false)
    _check_prefix(prefix)
    base = basename(prefix)
    ext = compress ? ".tsv.gz" : ".tsv"
    files = Dict{String,Any}()
    part(key, suffix) = (files[key] = base * suffix; prefix * suffix)
    write_bins(part("bins", "_bins.tsv"), m.grid)
    write_matrix(part("total", "_total" * ext), m)
    ntracks = m.allele === nothing ? 0 : length(m.allele)
    for h in 1:ntracks
        write_matrix(part("allele_$h", "_allele_$(h)" * ext), m; track = h)
    end
    g = m.grid
    meta = merge(_version_meta(), Dict{String,Any}(
        "format" => "CopyNumberPainter.matrix",
        "outputs" => sort!(collect(keys(files))), "files" => files,
        "assembly" => _assembly_meta(g.assembly),
        "grid" => Dict{String,Any}("size" => g.size, "max_masked_fraction" => g.max_masked_fraction,
                                   "nmasked" => g.nmasked,
                                   "mask_chrom" => [first(p) for p in g.mask],
                                   "mask_start" => [first(last(p)) for p in g.mask],
                                   "mask_stop" => [last(last(p)) for p in g.mask]),
        "rule" => string(nameof(typeof(m.rule))), "phasing" => String(m.phasing),
        "ncells" => ncells(m), "nbins" => nbins(g), "ntracks" => ntracks))
    open(io -> TOML.print(io, meta; sorted = true), prefix * "_meta.toml", "w")
    return prefix
end

"""
    load_matrix(prefix) -> CNMatrix

Read the files written by [`save_matrix`](@ref)`(prefix, m)`. The grid is rebuilt
from the metadata and must reproduce `<prefix>_bins.tsv` exactly.
"""
function load_matrix(prefix::AbstractString)
    meta = _check_bundle(TOML.parsefile(prefix * "_meta.toml"), "CopyNumberPainter.matrix", prefix)
    file(key) = joinpath(dirname(prefix), meta["files"][key])
    a = _assembly_from_meta(meta["assembly"])
    gm = meta["grid"]
    mask = Pair{String,UnitRange{Int}}[c => Int(s):Int(e) for (c, s, e) in
                                       zip(gm["mask_chrom"], gm["mask_start"], gm["mask_stop"])]
    g = BinGrid(a, Int(gm["size"]); mask = mask, max_masked_fraction = gm["max_masked_fraction"])
    read_bins(file("bins"), a) == g.bins || throw(ArgumentError(
        "$(prefix)_bins.tsv does not match the grid rebuilt from the metadata"))
    cells, names, total = _read_matrix_table(file("total"), nbins(g))
    ntracks = Int(meta["ntracks"])
    allele = ntracks == 0 ? nothing :
             [_read_matrix_table(file("allele_$h"), nbins(g))[3] for h in 1:ntracks]
    rule = get(_BIN_RULES, meta["rule"]) do
        throw(ArgumentError("unknown bin rule $(meta["rule"]) in $(prefix)_meta.toml"))
    end
    return CNMatrix(g, cells, names, total, allele, Symbol(meta["phasing"]), rule)
end
