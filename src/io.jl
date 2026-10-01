# Every writer takes either an IO or a path, and every reader an IO or a path; a path
# ending in `.gz` is (de)compressed on the fly. Tab-separated, one header line, `NA`
# for fields that do not apply. These are plain tables on purpose: the on-disk layout of
# a whole dataset directory is specified downstream, and inventing a second format here
# would guarantee a mismatch.
_with_io(f, io::IO) = f(io)
function _with_io(f, path::AbstractString)
    endswith(path, ".gz") || return open(f, path, "w")
    stream = CodecZlib.GzipCompressorStream(open(path, "w"))
    try
        return f(stream)
    finally
        close(stream)           # writes the gzip trailer; without it the file is truncated
    end
end

_with_input(f, io::IO) = f(io)
function _with_input(f, path::AbstractString)
    endswith(path, ".gz") || return open(f, path, "r")
    stream = CodecZlib.GzipDecompressorStream(open(path, "r"))
    try
        return f(stream)
    finally
        close(stream)
    end
end

# One tab-separated line, printed field by field. `join` built a String per row,
# which dominated the export time of million-row MEDICC2 files. The body is generated
# per argument-type tuple, so each field prints with a static call; a recursive
# version left the compiler to decide about inlining and measured twice as slow.
@generated function _printrow(io::IO, fields...)
    body = Expr(:block)
    for i in 1:length(fields)
        push!(body.args, :(print(io, fields[$i])))
        push!(body.args, :(write(io, $(i < length(fields) ? UInt8('\t') : UInt8('\n')))))
    end
    push!(body.args, :(nothing))
    return body
end

# Rows are collected in memory and handed on in ~1 MiB blocks. Many tiny writes are
# what makes a compressed stream, or a network filesystem, slow.
_flush_full!(io::IO, buf::IOBuffer) =
    (position(buf) >= 1 << 20 && write(io, take!(buf)); nothing)

# Column names, shared by the writers and the readers.
const _PROFILE_COLUMNS = ("node_id", "name", "chrom", "haplotype", "start", "stop", "cn")
const _EVENT_COLUMNS = ("node_id", "name", "order", "type", "chrom", "haplotype",
                        "start", "stop", "delta", "scale", "mode")
const _BIN_COLUMNS = ("bin_index", "chrom", "start", "stop")
const _TREE_COLUMNS = ("node_id", "parent", "label", "source_id", "birthtime",
                       "edge_divisions", "edge_mutations")

# A name printed into a TSV must not contain the separators that would split its row.
function _tsv_field(s::AbstractString)
    any(c -> c in ('\t', '\n', '\r'), s) && throw(ArgumentError(
        "name $(repr(s)) contains a tab or line break and cannot be written to a TSV table"))
    return String(s)
end

_na(x) = x === nothing ? "NA" : x
_parse_na(::Type{String}, s) = s == "NA" ? nothing : String(s)
_parse_na(::Type{T}, s) where {T} = s == "NA" ? nothing : parse(T, s)

# Rows of a table whose header must be exactly `cols`, split on tabs.
function _read_table(src, cols::Tuple)
    _with_input(src) do io
        lines = eachline(io)
        first_line = iterate(lines)
        first_line === nothing && throw(ArgumentError("empty table; expected header $(join(cols, '\t'))"))
        header = first(first_line)
        Tuple(split(header, '\t')) == cols || throw(ArgumentError(
            "unexpected header $(repr(header)); expected $(join(cols, '\t'))"))
        rows = Vector{Vector{SubString{String}}}()
        for (k, line) in enumerate(lines)
            isempty(line) && continue
            f = split(line, '\t')
            length(f) == length(cols) || throw(ArgumentError(
                "row $(k + 1) has $(length(f)) fields, expected $(length(cols))"))
            push!(rows, f)
        end
        return rows
    end
end

"""
    write_profiles(io_or_path, res; cells = 1:nnodes(res.tree))

Write copy-number segmentations as a long table with columns
`node_id, name, chrom, haplotype, start, stop, cn`.

Coordinates are this package's internal 1-based inclusive convention. `cells` selects
which nodes to write and defaults to **every** node, so the file contains the ancestral
truth as well as the observable tips. After a `retain_internal = false` run the missing
profiles are reconstructed with [`replay`](@ref), so the file is identical to that of a
full run.
"""
function write_profiles(dest, res::CNAEvolution;
                        cells::AbstractVector{<:Integer} = 1:nnodes(res.tree))
    return _write_profile_table(dest, res.tree, res.assembly, cells, _profiles_for(res, cells))
end

# The long profile table for `ids` and their `profs`; shared by write_profiles and the
# root file of a saved simulation.
function _write_profile_table(dest, tree::PhyloTree, a::GenomeAssembly, ids, profs)
    chrom = [chromname(a, c) for c in 1:nchromosomes(a)]
    # Names are checked before the destination is opened, so a refused name cannot
    # truncate an existing file.
    names = [_tsv_field(cellname(tree, id)) for id in ids]
    _with_io(dest) do io
        buf = IOBuffer()
        _printrow(buf, _PROFILE_COLUMNS...)
        for (id, p, nm) in zip(ids, profs, names)
            for c in 1:nchromosomes(a), h in 1:ploidy(a, c)
                for sg in p.segments[slot(a, c, h)]
                    _printrow(buf, id, nm, chrom[c], h, sg.start, sg.stop, sg.cn)
                end
            end
            _flush_full!(io, buf)
        end
        write(io, take!(buf))
    end
    return dest
end

"""
    write_events(io_or_path, res)

Write the complete event log as a table with columns
`node_id, name, order, type, chrom, haplotype, start, stop, delta, scale, mode`.

`type` is `segmental` or `wgd`. Columns that do not apply to a row hold `NA`:
a doubling has no chromosome, span, delta or scale, and a segmental alteration has no
mode. `node_id` identifies the edge by its child; the root's rows are truncal
alterations.
"""
function write_events(dest, res::CNAEvolution)
    a = res.assembly
    chrom = [chromname(a, c) for c in 1:nchromosomes(a)]
    names = [_tsv_field(cellname(res.tree, le.node)) for le in res.events]
    _with_io(dest) do io
        buf = IOBuffer()
        _printrow(buf, _EVENT_COLUMNS...)
        for (le, nm) in zip(res.events, names)
            e = le.event
            if e isa SegmentalCNA
                _printrow(buf, le.node, nm, le.order, "segmental", chrom[e.chrom],
                          e.haplotype, e.start, e.stop, e.delta, e.scale, "NA")
            else
                _printrow(buf, le.node, nm, le.order, "wgd", "NA", "NA",
                          "NA", "NA", "NA", "NA", e.mode)
            end
            _flush_full!(io, buf)
        end
        write(io, take!(buf))
    end
    return dest
end

"""
    write_bins(io_or_path, grid)

Write the bin manifest as `bin_index, chrom, start, stop`, in 1-based inclusive
coordinates. Pairs with a [`CNMatrix`](@ref) to say what its columns mean.
"""
function write_bins(dest, g::BinGrid)
    a = g.assembly
    chrom = [chromname(a, c) for c in 1:nchromosomes(a)]
    _with_io(dest) do io
        buf = IOBuffer()
        _printrow(buf, _BIN_COLUMNS...)
        for (i, b) in enumerate(g.bins)
            _printrow(buf, i, chrom[b.chrom], b.start, b.stop)
            _flush_full!(io, buf)
        end
        write(io, take!(buf))
    end
    return dest
end

"""
    write_medicc2(io_or_path, m; include_xy = false, normal_name = "diploid")

Write a [`CNMatrix`](@ref) as MEDICC2 input: a long TSV with columns
`sample_id, chrom, start, end, cn_a, cn_b`.

Conversions and conventions, all of them things MEDICC2 requires:

- **0-based half-open (BED) coordinates**, converted from this package's 1-based
  inclusive segments at this boundary and nowhere else.
- **Identical segmentation across every sample**, which holds by construction because
  every cell is projected onto the same grid.
- A reference sample named `normal_name` holding the assembly's starting karyotype
  (`1/1` on a diploid chromosome, `1/0` on a hemizygous one), which is the root MEDICC2
  measures distances from.
- **Autosomes only** by default, matching MEDICC2's own bulk analyses. Set
  `include_xy = true` to include the sex chromosomes; note a hemizygous chromosome
  exports as `cn_b = 0`.

Warns if any per-allele copy number (`cn_a` or `cn_b`) exceeds 8, which MEDICC2's
alphabet cannot represent; the total may exceed 8.

Only the observable cells belong in this file. The ancestral profiles are the ground
truth you compare MEDICC2's reconstruction *against* — write those with
[`write_profiles`](@ref) instead, and never feed them in as input.
"""
function write_medicc2(dest, m::CNMatrix; include_xy::Bool = false,
                       normal_name::AbstractString = "diploid")
    m.allele === nothing && throw(ArgumentError(
        "write_medicc2 needs allele-specific tracks; build the matrix with allele = true"))
    length(m.allele) == 2 || throw(ArgumentError(
        "MEDICC2 takes exactly two allele columns, but this matrix has " *
        "$(length(m.allele)) haplotype tracks; a third track would be silently dropped"))
    g = m.grid
    a = g.assembly
    keep = include_xy ? collect(1:nchromosomes(a)) : autosomes(a)
    cols = Int[]
    for c in keep
        append!(cols, bins_of(g, c))
    end
    A, B = m.allele[1], m.allele[2]
    # MEDICC2's limit is per allele (cn_a and cn_b each in 0..8), so a total above 8 is
    # fine as long as neither allele track exceeds it.
    allele_max = isempty(cols) || ncells(m) == 0 ? 0 : max(maximum(view(A, :, cols)), maximum(view(B, :, cols)))
    if allele_max > 8
        @warn "per-allele copy numbers above 8 cannot be represented by MEDICC2 (its alphabet caps at 8 per allele); the exported values will be out of range" max_allele_cn = allele_max
    end
    chrom = [chromname(a, c) for c in 1:nchromosomes(a)]
    normal = _tsv_field(normal_name)
    names = [_tsv_field(n) for n in m.names]
    # MEDICC2 would merge a cell sharing the normal's name with the normal sample.
    normal in names && throw(ArgumentError(
        "a cell is named $(repr(normal)), the same as the normal sample; pass a different normal_name"))
    _with_io(dest) do io
        buf = IOBuffer()
        _printrow(buf, "sample_id", "chrom", "start", "end", "cn_a", "cn_b")
        # The normal is the karyotype the assembly starts from, so a hemizygous
        # chromosome is 1/0 here too; a 1/1 normal would put a spurious loss on
        # every cell's branch.
        for col in cols
            b = g.bins[col]
            pl = ploidy(a, b.chrom)
            _printrow(buf, normal, chrom[b.chrom], b.start - 1, b.stop,
                      min(pl, 1), pl >= 2 ? 1 : 0)
        end
        _flush_full!(io, buf)
        for row in 1:ncells(m)
            nm = names[row]
            for col in cols
                b = g.bins[col]
                _printrow(buf, nm, chrom[b.chrom], b.start - 1, b.stop,
                          A[row, col], B[row, col])
            end
            _flush_full!(io, buf)
        end
        write(io, take!(buf))
    end
    return dest
end

# A node label as the tree table prints it; refuses what would not read back.
function _tree_label_field(t::PhyloTree, i::Integer)
    l = node(t, i).label
    l === nothing && return "NA"
    l == "NA" && throw(ArgumentError("node $i is labelled \"NA\", which reads back as a missing label"))
    return _tsv_field(l)
end

"""
    write_tree(io_or_path, tree)

Write `tree` as a node table with columns
`node_id, parent, label, source_id, birthtime, edge_divisions, edge_mutations`, rows
in preorder and `NA` where a field is `nothing`.

Unlike newick, which carries one number per edge, this is lossless: every field and
the order of every node's children survive [`read_tree`](@ref). A label spelled `NA`
cannot be told apart from a missing one and is refused.
"""
function write_tree(dest, t::PhyloTree)
    order = preorder(t)
    labels = [_tree_label_field(t, i) for i in order]
    _with_io(dest) do io
        buf = IOBuffer()
        _printrow(buf, _TREE_COLUMNS...)
        for (i, lb) in zip(order, labels)
            nd = node(t, i)
            _printrow(buf, i, _na(nd.parent), lb, _na(nd.source_id), _na(nd.birthtime),
                      _na(nd.edge_divisions), _na(nd.edge_mutations))
            _flush_full!(io, buf)
        end
        write(io, take!(buf))
    end
    return dest
end

"""
    read_tree(io_or_path) -> PhyloTree

Read a table written by [`write_tree`](@ref).
"""
function read_tree(src)
    rows = _read_table(src, _TREE_COLUMNS)
    n = length(rows)
    parent = Vector{Union{Int,Nothing}}(nothing, n)
    kids = [Int[] for _ in 1:n]
    bt = Vector{Union{Float64,Nothing}}(nothing, n)
    dv = Vector{Union{Int,Nothing}}(nothing, n)
    mu = Vector{Union{Int,Nothing}}(nothing, n)
    lb = Vector{Union{String,Nothing}}(nothing, n)
    sid = Vector{Union{Int,Nothing}}(nothing, n)
    seen = falses(n)
    for r in rows                                  # preorder: children arrive in order
        i = parse(Int, r[1])
        (1 <= i <= n && !seen[i]) || throw(ArgumentError("tree table node ids must be 1:$n, each once"))
        seen[i] = true
        parent[i] = _parse_na(Int, r[2])
        if parent[i] !== nothing
            1 <= parent[i] <= n || throw(ArgumentError("node $i has parent $(parent[i]), which is out of range"))
            push!(kids[parent[i]], i)
        end
        lb[i], sid[i] = _parse_na(String, r[3]), _parse_na(Int, r[4])
        bt[i], dv[i], mu[i] = _parse_na(Float64, r[5]), _parse_na(Int, r[6]), _parse_na(Int, r[7])
    end
    return PhyloTree([PhyloNode(i, parent[i], kids[i], bt[i], dv[i], mu[i], lb[i], sid[i]) for i in 1:n])
end

"""
    read_profiles(io_or_path, assembly) -> Dict{Int,CNProfile}

Read a table written by [`write_profiles`](@ref): node id => profile. Every profile is
checked with [`check_invariants`](@ref), so an incomplete or hand-edited table fails
here rather than downstream.
"""
function read_profiles(src, a::GenomeAssembly)
    segs = Dict{Int,Vector{Vector{Segment}}}()
    for r in _read_table(src, _PROFILE_COLUMNS)
        id = parse(Int, r[1])
        slots = get!(() -> [Segment[] for _ in 1:nslots(a)], segs, id)
        s = slot(a, chromindex(a, r[3]), parse(Int, r[4]))
        push!(slots[s], Segment(parse(Int, r[5]), parse(Int, r[6]), parse(Int, r[7])))
    end
    out = Dict{Int,CNProfile}()
    for (id, s) in segs
        p = CNProfile(a, s)
        try
            check_invariants(p)
        catch e
            e isa ErrorException || rethrow()
            throw(ArgumentError("profile of node $id is invalid: $(e.msg)"))
        end
        out[id] = p
    end
    return out
end

"""
    read_events(io_or_path, assembly) -> Vector{LoggedEvent}

Read a table written by [`write_events`](@ref), in file order.
"""
function read_events(src, a::GenomeAssembly)
    out = LoggedEvent[]
    for (k, r) in enumerate(_read_table(src, _EVENT_COLUMNS))
        ev = if r[4] == "wgd"
            WholeGenomeDoubling(Symbol(r[11]))
        elseif r[4] == "segmental"
            c, h, stop = chromindex(a, r[5]), parse(Int, r[6]), parse(Int, r[8])
            (h >= 1 && h <= ploidy(a, c)) || throw(ArgumentError(
                "event row $k: chromosome $(r[5]) has no haplotype $h"))
            stop <= chromlength(a, c) || throw(ArgumentError(
                "event row $k: stop $stop is beyond the length of $(r[5]) ($(chromlength(a, c)))"))
            SegmentalCNA(c, h, parse(Int, r[7]), stop, parse(Int, r[9]), Symbol(r[10]))
        else
            throw(ArgumentError("unknown event type $(repr(r[4]))"))
        end
        push!(out, LoggedEvent(parse(Int, r[1]), parse(Int, r[3]), ev))
    end
    return out
end

"""
    read_bins(io_or_path, assembly) -> Vector{Bin}

Read a bin manifest written by [`write_bins`](@ref).
"""
function read_bins(src, a::GenomeAssembly)
    out = Bin[]
    for (k, r) in enumerate(_read_table(src, _BIN_COLUMNS))
        c, start, stop = chromindex(a, r[2]), parse(Int, r[3]), parse(Int, r[4])
        (1 <= start <= stop <= chromlength(a, c)) || throw(ArgumentError(
            "bin row $k: $start-$stop is not within 1:$(chromlength(a, c)) of $(r[2])"))
        push!(out, Bin(c, start, stop))
    end
    return out
end

"""
    write_matrix(io_or_path, m; track = :total)

Write one track of a [`CNMatrix`](@ref) as a wide table: one row per cell with
columns `node_id, name, bin_1, …, bin_N`, where `bin_k` is row `k` of
[`write_bins`](@ref). `track` is `:total` or a haplotype index (`1`, `2`, …).

The shape most single-cell tools load directly. [`save_matrix`](@ref) writes every
track together with the metadata needed to read them back.
"""
function write_matrix(dest, m::CNMatrix; track::Union{Symbol,Integer} = :total)
    X = track === :total ? m.total :
        track isa Integer && m.allele !== nothing && 1 <= track <= length(m.allele) ? m.allele[track] :
        throw(ArgumentError("track must be :total or a haplotype index of this matrix, got $(repr(track))"))
    # Names are checked before the destination is opened, so a refused name cannot
    # truncate an existing file.
    names = [_tsv_field(nm) for nm in m.names]
    _with_io(dest) do io
        buf = IOBuffer()
        print(buf, "node_id\tname")
        for k in 1:nbins(m.grid)
            print(buf, "\tbin_", k)
        end
        print(buf, '\n')
        for row in 1:ncells(m)
            print(buf, m.cells[row], '\t', names[row])
            for col in 1:size(X, 2)
                print(buf, '\t', X[row, col])
            end
            print(buf, '\n')
            _flush_full!(io, buf)
        end
        write(io, take!(buf))
    end
    return dest
end

# Read a wide table written by write_matrix: (cells, names, values).
function _read_matrix_table(src, nbins::Int)
    _with_input(src) do io
        lines = eachline(io)
        first_line = iterate(lines)
        first_line === nothing && throw(ArgumentError("empty matrix table"))
        header = split(first(first_line), '\t')
        header == vcat(["node_id", "name"], ["bin_$k" for k in 1:nbins]) ||
            throw(ArgumentError("matrix table header does not match a $nbins-bin grid"))
        cells, names, rows = Int[], String[], Vector{Vector{Int}}()
        for line in lines
            isempty(line) && continue
            f = split(line, '\t')
            length(f) == nbins + 2 || throw(ArgumentError(
                "matrix row has $(length(f)) fields, expected $(nbins + 2)"))
            push!(cells, parse(Int, f[1])); push!(names, String(f[2]))
            push!(rows, [parse(Int, x) for x in @view f[3:end]])
        end
        X = isempty(rows) ? zeros(Int, 0, nbins) : permutedims(reduce(hcat, rows))
        return cells, names, X
    end
end
