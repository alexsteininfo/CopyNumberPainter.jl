# Output

## Running the simulation

```julia
res = simulate_cnas(tree, assembly, model; seed = 1)
res = simulate_cnas(tree, assembly, model)   # a seed is drawn and recorded in res.seed
res = simulate_cnas(tree, assembly, model; seed = 1, retain_internal = false)
res = simulate_cnas(tree, assembly, model; seed = 1, rng_mode = :per_node)
```

[`simulate_cnas`](@ref) descends depth-first and iteratively, so tree depth costs heap
rather than stack and a 10⁵-deep tree cannot overflow it. Memory is dominated by stored
profiles: every node's by default, the leaves' with `retain_internal = false`. In a lean
run, the working copies along the current root-to-node path are the only other profiles
alive.

### `rng_mode`

`:global` threads one stream through the whole traversal — the straightforward choice.

`:per_node` gives each edge its own stream, derived from the seed and the node's
`source_id`. An edge then draws *identically* no matter which other edges exist, so
simulating on a sampled tree yields **exactly** the alterations of simulating on the
full tree and subsetting — provided that:

1. every node carries a `source_id` that is the same in both trees (converted
   `NonMarkovEvolution.jl` trees do; a node without one falls back to its dense id,
   which sampling renumbers);
2. the sampled tree keeps unary nodes (pruning without collapsing);
3. the doubling policy commutes: `NoWGD`, `RateWGD` (drawn per edge under this mode),
   or `ScheduledWGD(...; by = :source_id, allow_missing = true)`. `ExactlyNWGD` never
   commutes and warns; `ScheduledWGD` by dense id commutes only if the ids happen to
   coincide;
4. any custom rate, target, extent or kind draws only from the `rng` it is passed.

It needs a seed, which is drawn if you give neither `seed` nor `rng`; an explicit `rng` is an error.

Seeded streams are version-stable. A seed sets Xoshiro's state through a fixed integer mixer, not through Julia's seeding hash, so the same seed gives the same random numbers on every Julia version. What can still change between versions is how a *sampler* turns those numbers into draws (for example, Distributions' `Poisson`), so a bit-identical rerun also needs the same package versions. `save_simulation` records them.

## What comes back

[`CNAEvolution`](@ref) holds the tree, the assembly, the model, the profiles, the
complete event log, the rejection tally and the seed, plus `rng_mode` (how random
streams were derived), `root_base` (the root's state before any root-logged event) and
`event_ranges` (the contiguous run of the log belonging to each edge, which [`events_on`](@ref) slices). When a simulation is loaded from
disk, the model is represented as a [`ModelRecord`](@ref), which captures the original
parameters and configuration in TOML-safe form.

```julia
profile(res, i)            # one node's profile; throws if it was not retained
leaf_profiles(res)         # the observable cells, in leaves(res.tree) order
res.profiles               # indexed by node id; entries may be nothing
```

Profiles are retained for **every** node by default — tips and internal nodes alike —
because the ancestral states are the ground truth an inference method's reconstruction
gets compared against. `retain_internal = false` keeps only the leaves for very large
trees.

### The event log is the primitive

```julia
res.events                 # Vector{LoggedEvent}, complete, always
nevents(res)
events_on(res, i)          # alterations on the edge into node i
events_below(res, i)       # everything strictly below node i
replay(res)                # reconstruct every profile from the log
```

A [`LoggedEvent`](@ref) records the edge (by its child), the event's `order` within
that edge's sequence, and the [`CNAEvent`](@ref) itself. Truncal alterations are logged
against the root. Events are stored in preorder of node and then by `order`, so
[`replay`](@ref) is a single forward pass.

The log is complete regardless of `retain_internal`, and the retained profiles are a
cache: `replay` reconstructs every profile from the root state plus the log. Nothing is
lost by running lean, and the equality of the two is a tested property.

### Rejections

```julia
res.rejections             # Dict{Symbol,Int}, keyed by the reason for the redraw
rejection_count(res)       # the total
```

A non-empty tally means the realised alteration distribution is conditioned on
viability. Report it.

The key `:no_effect` counts proposals that fell on absent DNA (copy number 0) and were
redrawn. Every logged event changed the genome, so `nevents` and per-edge event counts
are counts of real alterations; the `:no_effect` redraws condition the *target*
distribution on hitting existing material.

## Bin projection

Real low-coverage single-cell data is called in fixed bins, and every inference method
consumes a cells × bins integer matrix.

```julia
grid = BinGrid(assembly, 500_000)        # 500 kb, the DLP+ default
nbins(grid)
bins_of(grid, 1)                         # chr1's column range
```

Bins run consecutively within each chromosome, chromosomes appear in assembly order,
and chromosomes with zero ploidy contribute none — so a female grid has no `chrY`
columns. The final bin of each chromosome is short whenever the bin size does not
divide the chromosome length; that is documented rather than padded.

```julia
tot, alleles = project(profile(res, i), grid)
mat = CNMatrix(res, grid)                            # leaves by default
truth = CNMatrix(res, grid; cells = internal_nodes(res.tree))
mat = CNMatrix(res, grid; rule = AreaWeightedMean())
mat = CNMatrix(res, grid; allele = false)
ncells(mat); max_cn(mat)
```

[`CNMatrix`](@ref) carries `cells` (the node ids of its rows), `names` (from
[`cellname`](@ref), so rows match an exported newick's leaf labels), `total`, and
`allele` — one matrix per haplotype index — and `phasing`, which is `:haplotype` for
simulated (phased) data or `:major_minor` after [`major_minor`](@ref). A [`Bin`](@ref) on a chromosome whose ploidy
is below a haplotype index gets 0 on that track, which is how a male `chrX` exports
with `cn_b = 0`.

`total` is defined as the **sum of the haplotype tracks**, not as a reprojection of the
summed segmentation. Projection is non-linear, so the two differ at straddling bins;
defining it additively guarantees `total == A + B`, which every consumer of an
allele-specific matrix relies on.

### The straddling-bin rule

A bin containing a breakpoint has no unambiguously correct copy number. Two
[`BinRule`](@ref)s are provided:

- [`LengthWeightedMajority`](@ref) — the copy number covering the most base pairs in
  the bin, ties going to the lower value. The default.
- [`AreaWeightedMean`](@ref) — the length-weighted mean, rounded with halves up.

Either introduces a small systematic difference from a real caller's own binning.
Which one best matches a given caller is **unresolved** — see
[Limitations and open questions](limitations.md).

### Masking

Real callers filter bins in regions they cannot reliably call — centromeres, assembly
gaps, repetitive blacklists. The same observation model applies: regions can be masked
at the `BinGrid` level.

```julia
grid = BinGrid(assembly, 500_000; mask = centromere_mask(assembly))
grid = BinGrid(assembly, 500_000; mask = read_bed_mask("hg38-blacklist.bed"))
grid = BinGrid(assembly, 500_000; mask = vcat(centromere_mask(assembly), 
                                              read_bed_mask("hg38-blacklist.bed")))
```

A bin whose masked fraction exceeds `max_masked_fraction` (default 0.5, meaning 50%) is
dropped from the grid. Its column is then absent from every matrix and export, exactly
as if the caller had filtered it. The simulation itself is **unaffected**: masking
changes only what is observed, not the mutations themselves.

## Files

```julia
write_profiles("truth.tsv", res)                # node_id, name, chrom, haplotype, start, stop, cn
write_events("events.tsv", res)                 # node_id, name, order, type, chrom, haplotype, start, stop, delta, scale, mode
write_bins("bins.tsv", grid)                    # bin_index, chrom, start, stop
write_medicc2("cells.tsv", mat)                 # see the interoperability page
write_newick("tree.nwk", tree; branchlength = :divisions)
write_tree("tree.tsv", tree)                    # node_id, parent, label, source_id, birthtime, edge_divisions, edge_mutations
```

All writers accept an `IO` as well as a path. A path ending in `.gz` is written gzip-compressed, e.g.
`write_medicc2("cells.tsv.gz", mat)`; MEDICC2 reads it directly, since pandas infers
compression from the extension. Tables are tab-separated with one header
line and `NA` where a field does not apply — a doubling has no chromosome or span, a
segmental alteration has no mode. [`write_profiles`](@ref) writes every node by
default, so the file contains the ancestral truth as well as the tips;
[`write_events`](@ref) writes the whole log; [`write_bins`](@ref) says what a matrix's
columns mean.

These are deliberately plain tables. The on-disk layout of a whole dataset directory is
specified downstream, where the many-methods-many-datasets requirement lives, and
inventing a second format here would guarantee a mismatch. Newick is the interchange
format for trees.

Every writer has a reader. Readers take an `IO` or a path, and decompress a `.gz` path.

| Writer | Reader | Returns |
|:--|:--|:--|
| [`write_profiles`](@ref) | [`read_profiles`](@ref)`(src, assembly)` | `Dict{Int,CNProfile}`, node id to profile, validated |
| [`write_events`](@ref) | [`read_events`](@ref)`(src, assembly)` | `Vector{LoggedEvent}` in file order |
| [`write_bins`](@ref) | [`read_bins`](@ref)`(src, assembly)` | `Vector{Bin}` |
| [`write_tree`](@ref) | [`read_tree`](@ref)`(src)` | `PhyloTree` |
| [`write_newick`](@ref) | [`read_newick`](@ref) | `PhyloTree` |

[`write_tree`](@ref) is the lossless tree format: it keeps every node field and the
order of every node's children, which newick cannot, since newick carries a single
branch-length number. Use newick for interchange with other tools. Every writer refuses
a name containing a tab or line break, which would corrupt a row, and `write_tree`
refuses a label spelled `NA`, which would read back as missing.

## Saving and loading a simulation

The writers above produce one table each. To keep a whole run, [`save_simulation`](@ref)
writes a bundle of files sharing a prefix, and [`load_simulation`](@ref) reads it back:

```julia
save_simulation("out/run1", res)                 # out/run1_meta.toml, _tree.tsv, _events.tsv, _root.tsv, _profiles.tsv, .nwk
res2 = load_simulation("out/run1")               # CNAEvolution{ModelRecord}, complete
save_simulation("out/run1", res; profiles = :leaves, compress = true)
```

| File | Content |
|:--|:--|
| `<prefix>_meta.toml` | seed, `rng_mode`, versions, assembly, model description, rejection tally, counts, file list |
| `<prefix>_tree.tsv` | the lossless tree ([`write_tree`](@ref)) |
| `<prefix>_events.tsv[.gz]` | the complete event log |
| `<prefix>_root.tsv` | the root state before any root-logged event |
| `<prefix>_profiles.tsv[.gz]` | node profiles: every node (`profiles = :all`), the tips (`:leaves`) or none (`:none`) |
| `<prefix>.nwk` | the tree in newick, with the first branch-length field every edge has (`newick = :auto`) |

The metadata file records `format`, `format_version`, `package_version`,
`julia_version`, `outputs`, `files`, `assembly` and, for a simulation, `seed` (when
known), `rng_mode`, `retain_internal`, `nnodes`, `nleaves`, `nevents`, `saved_profiles`, `rejections`,
`newick_branchlength`, `model` and `model_repr`.

The event log is the primitive: a profile that was not saved, or whose table is missing
or unreadable (which gives a warning), is rebuilt by replaying the log, so a loaded
result is always complete and `retain_internal` is `true`. A model can hold closures
that a file cannot rebuild, so a loaded result carries a [`ModelRecord`](@ref), the
recorded description and printed form of what ran. Everything computed from the result
works as for a fresh one. `compress = true` gzips the event and profile tables.

A matrix has its own bundle. [`save_matrix`](@ref) writes `<prefix>_meta.toml` (assembly,
grid size and mask, bin rule, phasing), `<prefix>_bins.tsv`, `<prefix>_total.tsv[.gz]`
and one `<prefix>_allele_<h>.tsv[.gz]` per haplotype track; [`load_matrix`](@ref)
rebuilds the grid from the metadata and checks it against the bin manifest.
[`write_matrix`](@ref)`(dest, m; track = :total)` writes a single track as a wide table
`node_id, name, bin_1, …, bin_N`, the shape most single-cell tools load directly.

The TSV columns are the stable interface for tools outside Julia, and `format_version`
increases whenever a column changes. A bundle with a newer `format_version` than the
package reads is refused.
