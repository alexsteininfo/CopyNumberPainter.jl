# Using this as a MEDICC2 benchmark

This page collects what you need to simulate data that MEDICC2 can be scored on, and
what the ground truth does and does not tell you. The export format itself is on the
[MEDICC2 interoperability](interop.md) page.

## Settings that match MEDICC2's model

MEDICC2 counts a whole-genome doubling as *one event that adds one copy of every
allele still present*: copy number 0 is absorbing, so absent segments stay absent.
Later events then act on the doubled genome. Match that, or the
truth and the method are measuring different things:

```julia
using CopyNumberEvolution

# `tree` is the lineage tree built in the README quickstart.
assembly = hg38(:female)

model = CNAModel(
    rate    = PerDivision(0.5),
    wgd     = ScheduledWGD(mrca(tree, [4, 5]) => 1; mode = :increment),
    initial = TruncalCNAs(4; mode = :increment),   # the default `nothing` inherits the policy's mode
)

res  = simulate_cnas(tree, assembly, model; seed = 1)
grid = BinGrid(assembly, 500_000; mask = centromere_mask(assembly))
mat  = CNMatrix(res, grid)

maximum(maximum, mat.allele) <= 8 || @warn "MEDICC2's alphabet stops at 8 per allele"
```

- **`mode = :increment`** for every doubling, in the policy and in
  [`TruncalCNAs`](@ref). `WholeGenomeDoubling(:increment)` is MEDICC2's definition.
- **Each allele at most 8.** MEDICC2's alphabet caps `cn_a` and `cn_b` separately, so
  the total (`max_cn(mat)`) may exceed 8. Lower the rate or the number of doublings
  rather than capping in the simulation.
- **Autosomes.** Sex chromosomes need `include_xy` on export and a normal that is
  hemizygous where the genome is (see the interoperability page).
- **A 500 kb grid with `centromere_mask`**, which matches DLP+-style single-cell
  binning. Add [`read_bed_mask`](@ref) for a blacklist.

## Phased versus unphased input

MEDICC2 takes allele-specific profiles. Simulated data are phased, so two inputs are
possible from the same matrix:

```julia
write_medicc2("cells_phased.tsv.gz", mat)                    # haplotype-resolved
write_medicc2("cells_unphased.tsv.gz", major_minor(mat))     # major/minor per bin
```

[`write_medicc2`](@ref) writes a normal derived from the ploidy (1 allele on
hemizygous chromosomes, 0 on the absent one). The unphased input throws away which
haplotype an allele sits on, so a method is expected to do somewhat worse on it; the
gap is itself a useful benchmark.

## The ground truth to compare against

```julia
save_simulation("results/run", res)        # profiles of every node, the complete log, the tree
res = load_simulation("results/run")       # later, without re-simulating
events_on(res, i)                          # the events on the edge into node i (a view)
```

[`save_simulation`](@ref) keeps the profile of every node, including the internal ones
no method sees, and the complete event log. `compress = true` gzips the tables
(`save_simulation` takes a prefix, not a file name; `save_matrix` and `load_matrix` do
the same for a bin matrix).

Every logged event changed the genome: a proposal on DNA that is already absent is
redrawn and only tallied as `:no_effect`. The number of logged events on an edge is
therefore the true branch length in events, with no correction for events that did
nothing. It equals a MEDICC2-style event count only under `GainLoss(p, 1)` (a
copy-number change of exactly 1 per event) and `:increment` doublings; with larger
deltas one logged event corresponds to several unit changes.

## Known differences between truth and what any method can see

- **Events hidden by later ones.** A gain later overwritten by a loss of the same
  segment leaves no trace in the leaves, yet it is in the log. The distance between
  what the log counts and what the leaves support is the parsimony gap, which is what
  a minimum-evolution method is bound to underestimate.
- **Sub-bin focal events.** With the default length distribution, roughly a fifth to a
  quarter of focal events are shorter than a 500 kb bin. Projection onto the grid can blur or
  drop them, so they are in the truth and often absent from the input.
- **Unary chains versus MEDICC2's tree.** The simulated lineage tree has unary nodes
  (single-child divisions), which MEDICC2's output tree does not. To compare edge
  lengths, sum the true event counts along each unary chain. The diploid to MRCA edge
  is special: it includes the truncal events *and* everything that happened between
  the founder and the MRCA.
