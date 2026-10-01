# News

## 0.2.0

The package is renamed **CopyNumberEvolution.jl → CopyNumberPainter.jl**: replace
`using CopyNumberEvolution` with `using CopyNumberPainter`. The UUID is unchanged, so
existing environments switch by re-adding the package under its new name. Bundles
written by `save_simulation`/`save_matrix` under the old name still load.

Breaking changes to simulated output: the same seed gives different results than in 0.1.

- Focal events now start uniformly over every position from which they overlap the
  chromosome and are clipped at both ends, so coverage is uniform along it. In 0.1,
  q-telomeres were covered up to 6.6× more often than p-telomeres.
- Whole-genome doublings fall at a uniformly random position among their edge's
  alterations instead of always first. `TruncalCNAs(n; wgd, mode)` places truncal doublings.
- `CNWeighted(β = 1.0; length_weighted = true)` weights slots by copy-number material,
  length × mean copy number^β. `length_weighted = false` restores 0.1 behaviour.
- `read_newick(...; branchlength = :time)` no longer sets `edge_divisions = 1`.
- Proposals falling entirely on absent DNA are redrawn (tallied as `:no_effect`)
  instead of being logged as events that change nothing, so every logged event is a
  real alteration. `rejection_count` therefore includes these `:no_effect` redraws.
- `CNAModel()`'s default target is now `CNWeighted()`; it was `UniformChromosome()`.

Changed:

- The simulator bridge follows the upstream rename of MutationLoadDynamics.jl to
  NonMarkovEvolution.jl (compat `NonMarkovEvolution = "0.4"`): the package extension is now
  `CopyNumberPainterNonMarkovEvolutionExt`.
- Version-stable seeding: SplitMix64 expands a seed into Xoshiro state words, so the same
  seed gives the same streams on every Julia version. Results differ from 0.1.
- A seed is drawn and recorded when none is given.
- `rng_mode`, `root_base` and `event_ranges` are stored on the result, so `replay` needs
  no model.
- `events_on` returns a view.
- An ambiguous `node_by_label` and a duplicate `source_id` are errors.
- The event log is validated (`_event_ranges`); a malformed one is refused.
- MEDICC2's cap of 8 is checked per allele, not per total.
- `mean_cn` is exported; `CNMatrix` has a `rule` field.
- `write_newick(...; branchlength = :none)` omits branch lengths.

New:

- `BinGrid(...; mask, max_masked_fraction)`, `centromere_mask`, `read_bed_mask`.
- `major_minor(mat)` for unphased (major/minor) export; the `CNMatrix.phasing` field.
- `ScheduledWGD(...; by = :source_id, allow_missing)`; under `rng_mode = :per_node`,
  `RateWGD` is drawn per edge and commutes with sampling; `ExactlyNWGD` warns.
- `initial_profile(init, assembly)`: the whole interface of a custom `InitialState`.
- `save_simulation`/`load_simulation`: a TOML metadata file plus TSV parts (`.gz` with
  `compress = true`); missing profiles are rebuilt by replay. `save_matrix`/`load_matrix`
  and `write_matrix` do the same for a `CNMatrix`.
- `write_tree`/`read_tree` (lossless node table); `read_profiles`, `read_events`, `read_bins`.
- Every `.gz` path is gzip-compressed or -decompressed.
- `ModelRecord`: the recorded model description stored with a saved simulation.
- `NodeRef`, with the AbstractTrees extension; `PhyloTree(::LeafSample)`.

Performance (1024 x 6073 matrix; approximate ranges, measured on a loaded shared host):

- `write_medicc2` 5.7 s to roughly 3-4 s; `CNMatrix` 0.95 s to roughly 0.4-0.6 s.

Fixed:

- `read_newick` no longer overflows the stack on deep trees.
- The MEDICC2 normal is `1/0` on hemizygous chromosomes (it was `1/1`); more than two
  haplotype tracks is an error instead of a silent drop.
- Duplicate output names are an error in `CNMatrix` and `write_newick`.
- `write_profiles` works after `retain_internal = false`.
- Weighted sampling can no longer select a zero weight; `CustomRate`, `Given` and
  `CNAModel` validate their inputs; labels with tabs or newlines are quoted.
