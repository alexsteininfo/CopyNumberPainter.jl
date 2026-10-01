"""
    Bin(chrom, start, stop)

One fixed-width window of a [`BinGrid`](@ref), 1-based inclusive. The last bin of each
chromosome is shorter than the rest whenever the bin size does not divide the
chromosome length.
"""
struct Bin
    chrom::Int
    start::Int
    stop::Int
end

Base.length(b::Bin) = b.stop - b.start + 1

"""
    BinGrid(assembly, size = 500_000; mask = [], max_masked_fraction = 0.5)

A fixed genomic binning — the resolution at which real low-coverage single-cell data
is called, and the form every downstream inference method consumes.

Bins run consecutively within each chromosome and chromosomes appear in assembly
order. Chromosomes with zero ploidy contribute no bins, so a female grid has no `chrY`
columns. The default 500 kb matches DLP+.

`mask` lists regions real callers drop — centromeres, assembly gaps, blacklists — as
`"chr" => start:stop` pairs in 1-based inclusive coordinates; see
[`centromere_mask`](@ref) and [`read_bed_mask`](@ref). A bin whose masked fraction
exceeds `max_masked_fraction` is removed from the grid, so its column is absent from
every matrix and export, as in a caller's filtered output. The simulation itself is
unaffected: masking changes only what is observed.

Fields: `assembly`, `size`, `bins::Vector{Bin}`, `chromranges::Vector{UnitRange{Int}}`
giving each chromosome's column range, and for provenance `mask` (merged and clipped),
`max_masked_fraction` and `nmasked`, the number of bins removed.

Its fields are part of the public API and are covered by semantic versioning.
"""
struct BinGrid
    assembly::GenomeAssembly
    size::Int
    bins::Vector{Bin}
    chromranges::Vector{UnitRange{Int}}
    mask::Vector{Pair{String,UnitRange{Int}}}
    max_masked_fraction::Float64
    nmasked::Int

    function BinGrid(a::GenomeAssembly, size::Integer = 500_000;
                     mask = Pair{String,UnitRange{Int}}[],
                     max_masked_fraction::Real = 0.5)
        size >= 1 || throw(ArgumentError("bin size must be ≥ 1 bp, got $size"))
        0 <= max_masked_fraction <= 1 || throw(ArgumentError(
            "max_masked_fraction must lie in [0, 1], got $max_masked_fraction"))
        merged = _merge_mask(a, mask)
        bins = Bin[]
        ranges = Vector{UnitRange{Int}}(undef, nchromosomes(a))
        nmasked = 0
        for c in 1:nchromosomes(a)
            first_i = length(bins) + 1
            if ploidy(a, c) == 0
                ranges[c] = first_i:(first_i - 1)   # empty
                continue
            end
            L = chromlength(a, c)
            iv = merged[c]
            k = 1                                   # first interval that may still overlap
            s = 1
            while s <= L
                e = min(L, s + size - 1)
                while k <= length(iv) && last(iv[k]) < s
                    k += 1
                end
                masked = 0
                j = k
                while j <= length(iv) && first(iv[j]) <= e
                    masked += min(e, last(iv[j])) - max(s, first(iv[j])) + 1
                    j += 1
                end
                if masked / (e - s + 1) > max_masked_fraction
                    nmasked += 1
                else
                    push!(bins, Bin(c, s, e))
                end
                s = e + 1
            end
            ranges[c] = first_i:length(bins)
        end
        flat = Pair{String,UnitRange{Int}}[chromname(a, c) => r
                                           for c in 1:nchromosomes(a) for r in merged[c]]
        new(a, Int(size), bins, ranges, flat, Float64(max_masked_fraction), nmasked)
    end
end

# Resolve `mask` against `a`: one sorted list of disjoint, clipped intervals per
# chromosome. Chromosomes the assembly does not know (chrM, alt contigs — common in
# published blacklists) are skipped with a single warning naming them.
function _merge_mask(a::GenomeAssembly, mask)
    per = [UnitRange{Int}[] for _ in 1:nchromosomes(a)]
    unknown = Set{String}()
    for (name, r) in mask
        c = findfirst(k -> chromname(a, k) == name, 1:nchromosomes(a))
        if c === nothing
            push!(unknown, String(name))
            continue
        end
        lo, hi = max(1, first(r)), min(chromlength(a, c), last(r))
        lo <= hi && push!(per[c], lo:hi)
    end
    isempty(unknown) ||
        @warn "mask names chromosomes not in assembly $(a.name); they were ignored" chromosomes = sort!(collect(unknown))
    for c in eachindex(per)
        out = UnitRange{Int}[]
        for r in sort!(per[c]; by = first)
            if !isempty(out) && first(r) <= last(out[end]) + 1
                out[end] = first(out[end]):max(last(out[end]), last(r))
            else
                push!(out, r)
            end
        end
        per[c] = out
    end
    return per
end

"""
    centromere_mask(assembly) -> Vector{Pair{String,UnitRange{Int}}}

Every chromosome's centromeric interval, as a `mask` for [`BinGrid`](@ref).
"""
centromere_mask(a::GenomeAssembly) =
    Pair{String,UnitRange{Int}}[chromname(a, c) => centromere(a, c) for c in 1:nchromosomes(a)]

"""
    read_bed_mask(path_or_io) -> Vector{Pair{String,UnitRange{Int}}}

Read masked regions from a BED file (the first three columns), converting its 0-based
half-open intervals to this package's 1-based inclusive ones. Blank, `#`, `track` and
`browser` lines are skipped (a first field equal to `track` or `browser`, so a contig
named `trackA` is kept), and fields may be separated by tabs or spaces. Zero-length
records (`start == end`) are rejected. A path ending in `.gz` is decompressed. Pass the
result as `mask` to [`BinGrid`](@ref).
"""
function read_bed_mask(io::IO)
    out = Pair{String,UnitRange{Int}}[]
    for (lineno, line) in enumerate(eachline(io))
        s = strip(line)
        isempty(s) && continue
        f = split(s)
        # Whole-field match: a contig called `trackA` is a region, not a header.
        (startswith(s, '#') || f[1] == "track" || f[1] == "browser") && continue
        length(f) >= 3 || throw(ArgumentError("BED line $lineno has fewer than 3 fields: $(repr(line))"))
        start0, stop = tryparse(Int, f[2]), tryparse(Int, f[3])
        (start0 === nothing || stop === nothing) &&
            throw(ArgumentError("BED line $lineno: start and end must be integers, got $(repr(line))"))
        0 <= start0 < stop ||
            throw(ArgumentError("BED line $lineno: need 0 ≤ start < end, got $start0 and $stop"))
        push!(out, String(f[1]) => (start0 + 1):stop)
    end
    return out
end
read_bed_mask(path::AbstractString) = _with_input(read_bed_mask, path)

"""
    nbins(grid) -> Int

Total number of bins across all chromosomes.
"""
nbins(g::BinGrid) = length(g.bins)

"""
    bins_of(grid, chrom) -> UnitRange{Int}

Column indices belonging to chromosome `chrom`; empty when it has no slots.
"""
bins_of(g::BinGrid, c::Integer) = g.chromranges[c]

Base.show(io::IO, g::BinGrid) = print(io, "BinGrid(", g.assembly.name, ", ", g.size,
    " bp, ", nbins(g), " bins", g.nmasked > 0 ? ", $(g.nmasked) masked" : "", ")")

"""
    BinRule

How to assign an integer copy number to a bin that straddles a breakpoint.

There is no unambiguously correct answer, and either rule introduces a small
systematic difference from a real caller's own binning.

!!! note "Open question"
    Which rule best matches a given caller's behaviour is unresolved. The default is a
    documented placeholder, not a settled decision — see the manual's Limitations page.

Concrete rules: [`LengthWeightedMajority`](@ref), [`AreaWeightedMean`](@ref).
"""
abstract type BinRule end

"""
    LengthWeightedMajority()

Give the bin the copy number covering the most base pairs within it. Ties resolve to
the **lower** copy number, so the rule is deterministic. The default.
"""
struct LengthWeightedMajority <: BinRule end

"""
    AreaWeightedMean()

Give the bin the length-weighted mean copy number, rounded to the nearest integer with
halves going **up**.
"""
struct AreaWeightedMean <: BinRule end

"""
    project(segs, grid, chrom, rule) -> Vector{Int}

Project one slot's segmentation onto the bins of chromosome `chrom`, returning one
integer per bin of that chromosome.
"""
function project(segs::Vector{Segment}, g::BinGrid, c::Integer, rule::BinRule)
    cols = bins_of(g, c)
    return _project_into!(Vector{Int}(undef, length(cols)), segs, g, cols, rule)
end

# Values of the bins `cols` (all on one chromosome, ascending) from one slot's
# segmentation, in a single forward sweep: segments and bins are both sorted, so each
# is visited a bounded number of times instead of one binary search per bin. Masked
# bins leave gaps between bins, so the pointer advances by bin start, never by count.
function _project_into!(dest::AbstractVector{Int}, segs::Vector{Segment}, g::BinGrid,
                        cols::AbstractUnitRange{Int}, rule::BinRule)
    i = 1
    for (k, col) in enumerate(cols)
        b = g.bins[col]
        while segs[i].stop < b.start
            i += 1
        end
        dest[k] = _bin_value(segs, i, b, rule)
    end
    return dest
end

# The value of bin `b`, given the index `i` of the segment containing b.start.
function _bin_value(segs::Vector{Segment}, i::Int, b::Bin, ::LengthWeightedMajority)
    bestcn, bestlen = segs[i].cn, 0
    while i <= length(segs) && segs[i].start <= b.stop
        sg = segs[i]
        overlap = min(sg.stop, b.stop) - max(sg.start, b.start) + 1
        if overlap > bestlen || (overlap == bestlen && sg.cn < bestcn)
            bestcn, bestlen = sg.cn, overlap
        end
        i += 1
    end
    return bestcn
end

function _bin_value(segs::Vector{Segment}, i::Int, b::Bin, ::AreaWeightedMean)
    acc = 0
    while i <= length(segs) && segs[i].start <= b.stop
        sg = segs[i]
        acc += (min(sg.stop, b.stop) - max(sg.start, b.start) + 1) * sg.cn
        i += 1
    end
    return round(Int, acc / length(b), RoundNearestTiesUp)
end

"""
    project(profile, grid; rule = LengthWeightedMajority()) -> (total, alleles)

Project a whole profile onto `grid`.

Returns `total::Vector{Int}`, one value per bin, and `alleles::Vector{Vector{Int}}`,
one track per haplotype index. A bin on a chromosome whose ploidy is below a haplotype
index gets 0 on that track, which is how a male `chrX` gets `cn_b = 0`.

`total` is the **sum of the haplotype tracks**, not a reprojection of the summed
segmentation. Bin projection is non-linear, so the two differ at straddling bins;
defining it additively guarantees `total == A + B`, which every consumer of an
allele-specific matrix relies on.
"""
function project(p::CNProfile, g::BinGrid; rule::BinRule = LengthWeightedMajority())
    a = p.assembly
    same_assembly(a, g.assembly) || throw(ArgumentError(
        "profile is on $(a.name)/:$(a.sex) but the grid is on $(g.assembly.name)/:$(g.assembly.sex)"))
    maxploidy = maximum(a.ploidy)
    alleles = [zeros(Int, nbins(g)) for _ in 1:maxploidy]
    for c in 1:nchromosomes(a)
        cols = bins_of(g, c)
        isempty(cols) && continue
        for h in 1:ploidy(a, c)
            _project_into!(view(alleles[h], cols), p.segments[slot(a, c, h)], g, cols, rule)
        end
    end
    total = zeros(Int, nbins(g))
    for track in alleles
        total .+= track
    end
    return (total, alleles)
end

"""
    CNMatrix

Copy-number profiles projected onto a bin grid: the cells × bins integer matrix that
inference methods consume.

# Fields
- `grid::BinGrid`.
- `cells::Vector{Int}` — the tree node ids of the rows, in row order.
- `names::Vector{String}` — output names of the rows, from [`cellname`](@ref), so a
  row's identity matches its leaf label in an exported newick tree.
- `total::Matrix{Int}` — cells × bins total copy number.
- `allele::Union{Nothing,Vector{Matrix{Int}}}` — one cells × bins matrix per haplotype
  index, or `nothing`. `total` is always their sum.
- `phasing::Symbol` — `:haplotype` when track `h` is haplotype `h`, as simulated (phased by construction), or `:major_minor` after [`major_minor`](@ref).
- `rule::BinRule` — the straddling-bin rule the matrix was projected with.

Its fields are part of the public API and are covered by semantic versioning.
"""
struct CNMatrix
    grid::BinGrid
    cells::Vector{Int}
    names::Vector{String}
    total::Matrix{Int}
    allele::Union{Nothing,Vector{Matrix{Int}}}
    phasing::Symbol
    rule::BinRule
end

"""
    CNMatrix(res, grid; cells = leaves(res.tree), rule = LengthWeightedMajority(), allele = true)

Project the profiles of `cells` from a [`CNAEvolution`](@ref) onto `grid`.

`cells` defaults to the tree's leaves — the observable cells. Pass internal node ids
to build the ground-truth ancestral matrix that an inference method's reconstruction
can be compared against; keep that separate from what you hand the method as input.
"""
function CNMatrix(res::CNAEvolution, g::BinGrid;
                  cells::AbstractVector{<:Integer} = leaves(res.tree),
                  rule::BinRule = LengthWeightedMajority(),
                  allele::Bool = true)
    same_assembly(res.assembly, g.assembly) ||
        throw(ArgumentError(
            "result is on $(res.assembly.name)/:$(res.assembly.sex) but the grid is on " *
            "$(g.assembly.name)/:$(g.assembly.sex)"))
    ids = collect(Int, cells)
    n = length(ids)
    nb = nbins(g)
    a = res.assembly
    maxploidy = maximum(a.ploidy)
    # Filled bins × cells, so each cell writes one contiguous column (Julia is
    # column-major), then transposed once to the documented cells × bins layout.
    tracksT = [zeros(Int, nb, n) for _ in 1:maxploidy]
    for (j, p) in enumerate(_profiles_for(res, ids))
        for c in 1:nchromosomes(a)
            cols = bins_of(g, c)
            isempty(cols) && continue
            for h in 1:ploidy(a, c)
                _project_into!(view(tracksT[h], cols, j), p.segments[slot(a, c, h)], g, cols, rule)
            end
        end
    end
    total = permutedims(reduce(+, tracksT))
    tracks = allele ? [permutedims(T) for T in tracksT] : nothing
    names = [cellname(res.tree, id) for id in ids]
    _check_unique_names(names, ids)
    return CNMatrix(g, ids, names, total, tracks, :haplotype, rule)
end

"""
    ncells(m) -> Int

Number of rows in a [`CNMatrix`](@ref).
"""
ncells(m::CNMatrix) = length(m.cells)

"""
    max_cn(m) -> Int

Largest *total* copy number anywhere in the matrix. MEDICC2's limit of 8 applies per
allele, not to the total, so a total above 8 is not by itself a problem; check
`maximum(maximum, m.allele)` against 8 before export.
"""
max_cn(m::CNMatrix) = isempty(m.total) ? 0 : maximum(m.total)

Base.show(io::IO, m::CNMatrix) = print(io, "CNMatrix(", ncells(m), " cells × ",
    nbins(m.grid), " bins, max cn ", max_cn(m),
    m.phasing === :major_minor ? ", major/minor" : "", ")")
"""
    major_minor(m) -> CNMatrix

The unphased view of a phased matrix: in every bin of every cell, track 1 becomes the
**major** (larger) allele and track 2 the **minor** one, so `write_medicc2` then
writes `cn_a` = major and `cn_b` = minor.

Simulation output is always phased. This is the form real allele-specific callers
report without a phasing step, so use it to benchmark a method on realistic input,
e.g. MEDICC2's own evolutionary phasing. The assignment is made independently per
bin, so which haplotype is "major" may switch along a chromosome — exactly the
information unphased data lacks. `total` is unchanged, and `m` is not modified.
"""
function major_minor(m::CNMatrix)
    m.allele === nothing && throw(ArgumentError(
        "major_minor needs allele-specific tracks; build the matrix with allele = true"))
    length(m.allele) == 2 || throw(ArgumentError(
        "major_minor needs exactly two haplotype tracks, found $(length(m.allele))"))
    A, B = m.allele
    return CNMatrix(m.grid, copy(m.cells), copy(m.names), copy(m.total),
                    [max.(A, B), min.(A, B)], :major_minor, m.rule)
end

