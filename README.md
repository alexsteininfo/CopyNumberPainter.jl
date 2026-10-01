# CopyNumberEvolution.jl

[![Docs (dev)](https://img.shields.io/badge/docs-dev-blue.svg)](https://alexander-stein.github.io/CopyNumberEvolution.jl/dev/)

Forward simulation of somatic copy-number alterations along a cell-lineage tree.

Give it a tree — simulated by `NonMarkovEvolution.jl` or read from a newick file — and
it draws copy-number alterations along the edges from a diploid or given root state,
returning the allele-specific profile of **every** node plus a **complete log** of the
events that produced it, then projects those profiles onto a fixed bin grid.

This is an **observation model**: it does not infer trees, estimate parameters, or
compute distances between profiles.

## Install

```julia
using Pkg
Pkg.add(url = "https://github.com/alexander-stein/CopyNumberEvolution.jl")
```

## Quickstart

```julia
using CopyNumberEvolution

assembly = hg38(:female)

# A small lineage tree: root -> two divisions -> four leaves. Real trees come from
# NonMarkovEvolution.jl or a newick file:
#   tree = read_newick("lineage.nwk"; branchlength = :divisions)
tree = phylotree([nothing, 1, 1, 2, 2, 3, 3];
                 edge_divisions = [nothing, 1, 1, 1, 1, 1, 1])

model = CNAModel(
    rate    = PerDivision(0.5),
    target  = CNWeighted(1.0),
    extent  = ExtentMixture(p_chromosome = 0.05, p_arm = 0.15),
    kind    = GainLoss(0.6),
    wgd     = ScheduledWGD(mrca(tree, [4, 5]) => 1),
    initial = TruncalCNAs(4),
)

res  = simulate_cnas(tree, assembly, model; seed = 20260904)
grid = BinGrid(assembly, 500_000)
mat  = CNMatrix(res, grid)

mkpath("results")
save_simulation("results/run", res)      # profiles of every node plus the complete log
write_medicc2("cells.tsv.gz", mat)       # input for the reference method
res = load_simulation("results/run")     # reload the saved run later
```

To benchmark MEDICC2 with these simulations, see the manual page "Using this as a
MEDICC2 benchmark".

## What it models

Four event types, following MEDICC2's evolutionary model: segmental **gains** and
**losses** of arbitrary extent on a named haplotype, **loss of heterozygosity** as a
loss reaching zero, and **whole-genome doubling**. Whole-arm and whole-chromosome
events get their own probabilities, because no continuous length distribution produces
them at a realistic rate. Copy number 0 is **absorbing** — absent DNA is never
regained.

Alteration counts per edge come from an injectable rule: **per division**, **per unit
real time**, or **from the edge's mutation count**. That choice is the point of the
package — whether copy-number alterations accrue per division or per unit time is the
Markov-versus-non-Markov question transposed from point mutations, and the input tree
carries both quantities so both are computable.

Whole-genome doublings can be pinned to an **exact edge**
(`ScheduledWGD(mrca(tree, leaves) => 1)`), fixed in number with random placement, or
drawn from a rate.

## Scope

| | |
|:---|:---|
| Here | tree types, copy-number profiles, the alteration process, bin projection, MEDICC2 export |
| `NonMarkovEvolution.jl` | lineage-tree simulation and leaf sampling (a **weak** dependency) |
| Downstream | tree inference, distances, estimators — these must never depend on a simulator |
| Study repositories | parameter grids, file naming, figures |

Out of scope for v1, deliberately: a karyotype backend for general structural variants
(kept behind `apply!` so it can be added later), copy-number-neutral events,
breakage–fusion–bridge cycles, read-depth noise.

## Documentation

Build locally:

```bash
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate(); include("docs/make.jl")'
```

The manual covers the concepts, the input layer, every model component and the
scientific consequences of its defaults, the output and projection rules, MEDICC2
interoperability, and the open questions. `literature/MEDICC2.md` holds detailed notes
on the reference method.

## Examples

Five runnable scripts in `examples/`, each self-contained:

```bash
julia --project examples/01_quickstart.jl               # tree -> profiles -> files
julia --project examples/02_per_division_vs_per_time.jl # the Markov-vs-not question
julia --project examples/03_subclonal_wgd.jl             # mrca-scheduled WGD, :multiply vs :increment
julia --project examples/04_viability_rejection.jl       # AllowAll vs RejectAndRedraw
julia --project examples/05_newick_and_medicc2.jl        # newick in, MEDICC2 export out
```

## Open questions

Documented rather than silently settled — see the manual's Limitations page:

- The **straddling-bin projection rule** is unresolved; the default is a placeholder.
- The **copy-number ceiling** (MEDICC2 caps at 8 per allele) is warned about on export but not
  capped in simulation.
- `RejectAndRedraw`'s default `min_total_cn = 1` forbids homozygous deletions of *any*
  size, not just whole-chromosome nullisomy.
- **Rejection sampling conditions the model**: the realised alteration distribution is
  not the proposal distribution whenever rejections occur.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

The package-extension tests skip themselves unless `NonMarkovEvolution.jl` is
loadable. It is only a weak dependency, so run them from the repository root in a
temporary environment (never `Pkg.develop` into the package's own project, which would
write it into `[deps]`):

```bash
julia -e 'using Pkg; Pkg.activate(temp = true); Pkg.develop(path = pwd());
          Pkg.develop(path = "../NonMarkovEvolution.jl");
          Pkg.add(["Test", "Distributions", "Aqua", "AbstractTrees", "REPL"]);
          include(joinpath(pwd(), "test", "runtests.jl"))'
```

**No real patient data belongs in this repository.** Fixtures are synthetic and small;
`data/` is in `.gitignore`.

## License

See [LICENSE](LICENSE).
