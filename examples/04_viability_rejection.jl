# Example 4: viability constraints, and what "rejection sampling conditions the
# model" means in practice.
#
# An alteration can drive a region, a whole chromosome, or the single X of a male
# karyotype to copy number 0. Real data contains no cells with whole-chromosome
# nullisomy, so an unconstrained process can generate profiles that could not exist.
# See docs/src/model.md, "Viability", and docs/src/limitations.md.
#
# Run with: julia --project examples/04_viability_rejection.jl

using CopyNumberEvolution

# A hemizygous assembly: chr1 is present in a single copy, so the only whole-chromosome
# loss available drives it straight to copy number 0 — the sharpest possible test case.
assembly = GenomeAssembly("hemi", :male,
    [ChromosomeSpec("chr1", 1_000_000, 400_001:600_000),
     ChromosomeSpec("chr2", 1_000_000, 400_001:600_000)],
    [1, 2])   # chr1 hemizygous, chr2 the usual two copies

tree = phylotree([nothing, 1, 1, 1]; edge_divisions = [nothing, 1, 1, 1])

# A high rate, generous chromosome-loss probability, and a UniformChromosome target
# so nullisomy is likely unless something prevents it. Uniform targeting keeps
# proposing losses on chr1 until it is lost, which is what makes the contrast
# with RejectAndRedraw visible.
settings = (rate = PerDivision(3.0), extent = ExtentMixture(p_chromosome = 0.5),
            kind = GainLoss(0.3), target = UniformChromosome())

nullisomic(res) = count(leaves(res.tree)) do l
    any(sg -> sg.cn == 0, slot_segments(profile(res, l), 1, 1))
end

res_allow  = simulate_cnas(tree, assembly, CNAModel(; settings..., viability = AllowAll());
                            seed = 1)
res_reject = simulate_cnas(tree, assembly, CNAModel(; settings..., viability = RejectAndRedraw());
                            seed = 1)

println("AllowAll():        ", nullisomic(res_allow), "/", length(leaves(tree)),
        " leaves have chr1 nullisomy; rejections = ", res_allow.rejections)
println("RejectAndRedraw():  ", nullisomic(res_reject), "/", length(leaves(tree)),
        " leaves have chr1 nullisomy; rejections = ", res_reject.rejections,
        " (", rejection_count(res_reject), " total)")

println("\nSame tree, same model, same seed — the only difference is the viability rule.",
        "\nRejectAndRedraw() blocks every proposal that would reduce any chromosome to total",
        "\ncopy number 0 (the `min_total_cn = 1` constraint) and redraws instead, so no",
        "\nimpossible cell reaches the output. But this makes the realised alteration",
        "\ndistribution conditional on viability, not the raw proposal distribution —",
        "\nreport `rejection_count(res)` whenever it is non-zero.",
        "\n\nNote that `min_total_cn = 1` forbids homozygous deletions of ANY size on a",
        "\ndiploid chromosome, not only whole-chromosome loss; see docs/src/limitations.md",
        "\nif focal biallelic loss should be permitted instead.")

# Next: example 05 reads a tree from a newick file and exports for MEDICC2.
