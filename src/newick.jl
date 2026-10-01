const BRANCHLENGTH_MODES = (:divisions, :mutations, :time)

# Mutable staging node used only during parsing.
mutable struct _RawNode
    parent::Union{Int,Nothing}
    children::Vector{Int}
    label::Union{String,Nothing}
    brlen::Union{Float64,Nothing}
end

"""
    parse_newick(str; branchlength) -> PhyloTree

Parse a newick string.

A newick file carries **one** number per edge, but this package needs both real time
and a division count, so you must say which the numbers are. There is deliberately no
default: an implicit convention is the kind of thing that becomes invisible and wrong.

| `branchlength` | field populated | rate rules it enables |
|:---|:---|:---|
| `:divisions` | `edge_divisions` | [`PerDivision`](@ref) |
| `:mutations` | `edge_mutations` | [`FromEdgeMutations`](@ref) |
| `:time` | `birthtime`, by cumulative sum from the root | [`PerTime`](@ref) |

A time tree carries no division count: an edge's elapsed time — the interdivision time on a lineage tree — is [`edge_time`](@ref), and `PerDivision` throws on it rather than assuming one division per edge.

Fields not implied by the chosen mode stay `nothing`, and each rate rule throws a
named error when the field it needs is absent.

Supported syntax: named and unnamed internal nodes, arbitrary arity, single-quoted
labels (with `''` for a literal quote), `[...]` comments, missing branch lengths, and
arbitrary whitespace. Under `:divisions` and `:mutations` a fractional branch length
is rounded and warned about, since a fractional count is usually a sign the file
should have been read as `:time`.

# Examples
```jldoctest
julia> t = parse_newick("((A:1.0,B:2.0)X:0.5,C:3.0)R;"; branchlength = :time);

julia> node(t, node_by_label(t, "A")).birthtime
1.5
```
"""
function parse_newick(str::AbstractString; branchlength::Symbol = :unset)
    branchlength in BRANCHLENGTH_MODES || throw(ArgumentError(
        "branchlength must be one of $BRANCHLENGTH_MODES; got :$branchlength. " *
        "There is no default on purpose — say whether the numbers are divisions, mutations or time."))
    raw = _parse_raw(str)
    return _raw_to_tree(raw, branchlength)
end

"""
    read_newick(path_or_io; branchlength) -> PhyloTree

Read a newick tree from a file path or an `IO`. See [`parse_newick`](@ref) for the
meaning of `branchlength`.
"""
read_newick(io::IO; branchlength::Symbol = :unset) =
    parse_newick(read(io, String); branchlength = branchlength)
read_newick(path::AbstractString; branchlength::Symbol = :unset) =
    _with_input(io -> read_newick(io; branchlength = branchlength), path)

const _NEWICK_DELIMS = ('(', ')', ',', ':', ';', '[')

# Read position over the input. A struct rather than closures over a local `pos`:
# a captured, reassigned variable is boxed by Julia, which makes the parser
# type-unstable.
mutable struct _Cursor
    s::Vector{Char}
    pos::Int
end

_eof(c::_Cursor) = c.pos > length(c.s)
_peek(c::_Cursor) = c.s[c.pos]

function _skipspace!(c::_Cursor)
    while !_eof(c)
        ch = _peek(c)
        if isspace(ch)
            c.pos += 1
        elseif ch == '['
            depth = 0
            while !_eof(c)
                _peek(c) == '[' && (depth += 1)
                _peek(c) == ']' && (depth -= 1)
                c.pos += 1
                depth == 0 && break
            end
            depth == 0 || throw(ArgumentError("unterminated comment in newick input"))
        else
            return
        end
    end
end

function _readlabel!(c::_Cursor)
    _skipspace!(c)
    _eof(c) && return nothing
    if _peek(c) == '\''
        c.pos += 1
        buf = Char[]
        while true
            _eof(c) && throw(ArgumentError("unterminated quoted label in newick input"))
            ch = _peek(c)
            if ch == '\''
                if c.pos + 1 <= length(c.s) && c.s[c.pos + 1] == '\''
                    push!(buf, '\'')
                    c.pos += 2
                else
                    c.pos += 1
                    break
                end
            else
                push!(buf, ch)
                c.pos += 1
            end
        end
        return String(buf)
    end
    start = c.pos
    while !_eof(c) && !(_peek(c) in _NEWICK_DELIMS) && !isspace(_peek(c))
        c.pos += 1
    end
    return c.pos > start ? String(c.s[start:c.pos - 1]) : nothing
end

function _readlength!(c::_Cursor)
    _skipspace!(c)
    (!_eof(c) && _peek(c) == ':') || return nothing
    c.pos += 1
    _skipspace!(c)
    start = c.pos
    while !_eof(c) && (isdigit(_peek(c)) || _peek(c) in ('.', '-', '+', 'e', 'E'))
        c.pos += 1
    end
    c.pos > start || throw(ArgumentError("expected a number after ':' in newick input"))
    txt = String(c.s[start:c.pos - 1])
    val = tryparse(Float64, txt)
    val === nothing && throw(ArgumentError("could not parse branch length '$txt'"))
    val >= 0 || throw(ArgumentError("negative branch length $val is not allowed"))
    return val
end

# A node's trailing label and branch length, read once its subtree is complete.
function _finish_node!(c::_Cursor, nodes::Vector{_RawNode}, i::Int)
    nodes[i].label = _readlabel!(c)
    nodes[i].brlen = _readlength!(c)
    return nothing
end

# Called right after '(' or ',': the next thing must begin a subtree.
function _check_element_start(c::_Cursor, after_open::Bool)
    _eof(c) && return nothing              # reported as an unexpected end on the next read
    ch = _peek(c)
    after_open && ch == ')' && throw(ArgumentError("empty branch set '()' in newick input"))
    ch in (',', ')') && throw(ArgumentError(
        "empty element in a newick branch set at position $(c.pos): a stray, leading or trailing comma"))
    return nothing
end

# Iterative descent. `open` holds the internal nodes whose ')' is still ahead, so
# nesting depth costs heap, not stack: a 10^5-deep caterpillar parses like any other.
function _parse_raw(str::AbstractString)
    c = _Cursor(collect(str), 1)
    nodes = _RawNode[]
    _skipspace!(c)
    _eof(c) && throw(ArgumentError("empty newick input"))
    _peek(c) == ';' && throw(ArgumentError("empty newick tree ';' has no root node"))

    open = Int[]
    parent = nothing                       # parent of the next subtree to begin
    done = false
    while !done
        _skipspace!(c)
        _eof(c) && throw(ArgumentError("unexpected end of newick input"))
        push!(nodes, _RawNode(parent, Int[], nothing, nothing))
        me = length(nodes)
        parent === nothing || push!(nodes[parent].children, me)
        if _peek(c) == '('
            c.pos += 1
            _skipspace!(c)
            _check_element_start(c, true)
            push!(open, me)
            parent = me
            continue                       # descend into the first child
        end
        _finish_node!(c, nodes, me)        # a leaf
        # Climb: each ')' completes an open node; a ',' starts its next child.
        while true
            if isempty(open)
                done = true
                break
            end
            _skipspace!(c)
            _eof(c) && throw(ArgumentError("unbalanced parentheses in newick input"))
            ch = _peek(c)
            if ch == ','
                c.pos += 1
                _skipspace!(c)
                _check_element_start(c, false)
                parent = open[end]
                break
            elseif ch == ')'
                c.pos += 1
                _finish_node!(c, nodes, pop!(open))
            else
                throw(ArgumentError("expected ',' or ')' in newick input, found '$ch'"))
            end
        end
    end

    _skipspace!(c)
    (!_eof(c) && _peek(c) == ';') ||
        throw(ArgumentError("newick input must end with ';'" *
            (_eof(c) ? "" : " but continues with '$(_peek(c))'")))
    c.pos += 1
    _skipspace!(c)
    _eof(c) || throw(ArgumentError("trailing content after ';' in newick input"))
    return nodes
end

function _raw_to_tree(raw::Vector{_RawNode}, mode::Symbol)
    n = length(raw)
    order = _raw_preorder(raw)
    birthtimes = Vector{Union{Float64,Nothing}}(nothing, n)
    divisions = Vector{Union{Int,Nothing}}(nothing, n)
    mutations = Vector{Union{Int,Nothing}}(nothing, n)
    fractional = 0

    for i in order
        nd = raw[i]
        if mode === :time
            if nd.parent === nothing
                birthtimes[i] = 0.0
            else
                pbt = birthtimes[nd.parent]
                birthtimes[i] = (pbt === nothing || nd.brlen === nothing) ? nothing : pbt + nd.brlen
            end
        elseif nd.parent !== nothing && nd.brlen !== nothing
            v = nd.brlen
            abs(v - round(v)) > 1e-9 && (fractional += 1)
            iv = round(Int, v)
            mode === :divisions ? (divisions[i] = iv) : (mutations[i] = iv)
        end
    end

    fractional > 0 && @warn "newick branch lengths under :$mode should be whole numbers; $fractional fractional value(s) were rounded. If these are real times, read the file with branchlength = :time."

    nodes = [PhyloNode(i, raw[i].parent, raw[i].children,
                       birthtimes[i], divisions[i], mutations[i], raw[i].label, nothing)
             for i in 1:n]
    return PhyloTree(nodes)
end

function _raw_preorder(raw::Vector{_RawNode})
    r = findfirst(nd -> nd.parent === nothing, raw)
    out = Int[]
    stack = [r]
    while !isempty(stack)
        i = pop!(stack)
        push!(out, i)
        append!(stack, raw[i].children)
    end
    return out
end

const _NEWICK_SPECIAL = ('(', ')', ',', ':', ';', '[', ']', '\'')

_quote_label(name::AbstractString) =
    any(c -> c in _NEWICK_SPECIAL || isspace(c), name) ? "'" * replace(name, "'" => "''") * "'" : name

function _emit_name(t::PhyloTree, i::Integer, labels::Symbol)
    labels === :label && return cellname(t, i)
    labels === :id && return "cell_$(i)"
    labels === :source_id && begin
        sid = node(t, i).source_id
        sid === nothing && throw(ArgumentError("node $i has no source_id, so labels = :source_id cannot name it"))
        return "cell_$(sid)"
    end
    throw(ArgumentError("labels must be :label, :id or :source_id; got :$labels"))
end

function _emit_length(t::PhyloTree, i::Integer, mode::Symbol)
    if mode === :time
        return edge_time(t, i)
    elseif mode === :divisions
        v = node(t, i).edge_divisions
        v === nothing && throw(ArgumentError(
            "cannot write branchlength = :divisions: node $i has no edge_divisions"))
        return v
    else
        v = node(t, i).edge_mutations
        v === nothing && throw(ArgumentError(
            "cannot write branchlength = :mutations: node $i has no edge_mutations"))
        return v
    end
end

"""
    newick_string(tree; branchlength, labels = :label) -> String

Render `tree` as a newick string.

`branchlength` chooses **which** quantity the single branch-length field carries —
`:time`, `:divisions` or `:mutations` — and throws if any non-root node lacks it.
`:none` writes the topology and names without lengths. The
choice is explicit on write for the same reason it is on read: newick cannot carry
both real time and a division count, and a silent convention is a bug waiting to
happen.

`labels` picks the emitted names: `:label` uses [`cellname`](@ref) (a node's label, or
`cell_<id>`), `:id` always uses `cell_<dense id>`, and `:source_id` uses
`cell_<upstream id>`.
"""
function newick_string(t::PhyloTree; branchlength::Symbol = :unset, labels::Symbol = :label)
    branchlength in (BRANCHLENGTH_MODES..., :none) || throw(ArgumentError(
        "branchlength must be one of $BRANCHLENGTH_MODES, or :none to write no lengths; got :$branchlength"))
    labels in (:label, :id, :source_id) || throw(ArgumentError(
        "labels must be :label, :id or :source_id; got :$labels"))
    # Leaves only: internal labels are often support values ("100") and legitimately repeat.
    lv = leaves(t)
    _check_unique_names([_emit_name(t, i, labels) for i in lv], lv)
    io = IOBuffer()
    _write_subtree(io, t, treeroot(t), branchlength, labels)
    print(io, ';')
    return String(take!(io))
end

# Iterative emission: a frame is (node, child index, already-opened flag), so a
# 10^5-deep tree writes without recursion.
function _write_subtree(io::IO, t::PhyloTree, start::Int, mode::Symbol, labels::Symbol)
    stack = Tuple{Int,Int}[(start, 0)]
    while !isempty(stack)
        i, k = pop!(stack)
        kids = childrenof(t, i)
        if k == 0 && !isempty(kids)
            print(io, '(')
            push!(stack, (i, 1))
            push!(stack, (kids[1], 0))
        elseif k > 0 && k < length(kids)
            print(io, ',')
            push!(stack, (i, k + 1))
            push!(stack, (kids[k + 1], 0))
        else
            !isempty(kids) && print(io, ')')
            print(io, _quote_label(_emit_name(t, i, labels)))
            if mode !== :none && !isroot(t, i)
                print(io, ':', _emit_length(t, i, mode))
            end
        end
    end
    return io
end

"""
    write_newick(io_or_path, tree; branchlength, labels = :label)

Write `tree` in newick format to an `IO` or a file path. See
[`newick_string`](@ref).
"""
write_newick(io::IO, t::PhyloTree; branchlength::Symbol = :unset, labels::Symbol = :label) =
    print(io, newick_string(t; branchlength = branchlength, labels = labels))

function write_newick(path::AbstractString, t::PhyloTree; branchlength::Symbol = :unset,
                      labels::Symbol = :label)
    # Rendered before the file is opened, so a refused name cannot truncate it.
    str = newick_string(t; branchlength = branchlength, labels = labels)
    _with_io(io -> print(io, str), path)
end
