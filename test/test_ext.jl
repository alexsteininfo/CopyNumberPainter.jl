@testset "package extension" begin
    @testset "the extension is wired up in Project.toml" begin
        proj = read(joinpath(@__DIR__, "..", "Project.toml"), String)
        @test occursin("[weakdeps]", proj)
        @test occursin("NonMarkovEvolution = \"7b855ee6-6887-412f-a571-26d20a5a92d7\"", proj)
        @test occursin("CopyNumberPainterNonMarkovEvolutionExt = \"NonMarkovEvolution\"", proj)
        # and never a hard dependency: the simulator must not be reachable from an
        # inference-only install
        deps = match(r"\[deps\](.*?)\n\["s, proj)
        @test deps !== nothing
        @test !occursin("NonMarkovEvolution", deps.captures[1])
        @test isfile(joinpath(@__DIR__, "..", "ext",
                              "CopyNumberPainterNonMarkovEvolutionExt.jl"))
    end

    nme_loaded = try
        @eval using NonMarkovEvolution
        true
    catch
        false
    end

    if !nme_loaded
        @info """NonMarkovEvolution.jl is not available, so the conversion tests are skipped.
                 Run them from the repository root in a temporary environment:
                 julia -e 'using Pkg; Pkg.activate(temp = true); Pkg.develop(path = pwd()); Pkg.develop(path = "../NonMarkovEvolution.jl"); Pkg.add(["Test", "Distributions", "Aqua", "AbstractTrees", "REPL"]); include(joinpath(pwd(), "test", "runtests.jl"))'"""
    else
        @testset "converting a lineage tree" begin
            # Build a small BinaryNode tree by hand: founder -> two daughters,
            # the left one dividing again.
            root = NonMarkovEvolution.BinaryNode(
                NonMarkovEvolution.NonMarkovCell(1, 0.0, 3, 3, 1.0))
            NonMarkovEvolution.left_child!(root,
                NonMarkovEvolution.NonMarkovCell(2, 1.5, 4, 7, 1.0))
            NonMarkovEvolution.right_child!(root,
                NonMarkovEvolution.NonMarkovCell(3, 1.5, 1, 4, 1.0))
            NonMarkovEvolution.left_child!(root.left,
                NonMarkovEvolution.NonMarkovCell(4, 2.25, 2, 9, 1.0))
            NonMarkovEvolution.right_child!(root.left,
                NonMarkovEvolution.NonMarkovCell(5, 2.75, 0, 7, 1.0))

            t = PhyloTree(root)
            @test nnodes(t) == 5
            @test isroot(t, treeroot(t))
            @test node(t, treeroot(t)).source_id == 1
            @test node(t, treeroot(t)).birthtime == 0.0
            @test node(t, treeroot(t)).edge_divisions === nothing
            @test node(t, treeroot(t)).edge_mutations === nothing

            # every non-root edge is exactly one division
            @test all(node(t, i).edge_divisions == 1 for i in 1:nnodes(t) if !isroot(t, i))

            i2 = node_by_source_id(t, 2)
            @test node(t, i2).edge_mutations == 4
            @test edge_time(t, i2) ≈ 1.5
            @test node(t, node_by_source_id(t, 5)).edge_mutations == 0
            @test edge_time(t, node_by_source_id(t, 5)) ≈ 1.25

            @test Set(node(t, l).source_id for l in leaves(t)) == Set([3, 4, 5])
            @test founder_mutations(root) == 3
        end

        @testset "all three rate rules work on a converted tree" begin
            root = NonMarkovEvolution.BinaryNode(
                NonMarkovEvolution.NonMarkovCell(1, 0.0, 2, 2, 1.0))
            NonMarkovEvolution.left_child!(root,
                NonMarkovEvolution.NonMarkovCell(2, 1.0, 5, 7, 1.0))
            NonMarkovEvolution.right_child!(root,
                NonMarkovEvolution.NonMarkovCell(3, 2.0, 7, 9, 1.0))
            t = PhyloTree(root)
            rng = Random.Xoshiro(1)
            i2 = node_by_source_id(t, 2)
            @test n_cnas(PerDivision(1.0), t, i2, rng) isa Int
            @test n_cnas(PerTime(1.0), t, i2, rng) isa Int
            @test n_cnas(FromEdgeMutations(), t, i2, rng) == 5
        end

        @testset "a unary chain, as pruning leaves behind, converts unchanged" begin
            root = NonMarkovEvolution.BinaryNode(
                NonMarkovEvolution.NonMarkovCell(1, 0.0, 0, 0, 1.0))
            NonMarkovEvolution.left_child!(root,
                NonMarkovEvolution.NonMarkovCell(2, 1.0, 1, 1, 1.0))
            NonMarkovEvolution.left_child!(root.left,
                NonMarkovEvolution.NonMarkovCell(3, 2.0, 1, 2, 1.0))
            t = PhyloTree(root)
            @test nnodes(t) == 3
            @test leaves(t) == [node_by_source_id(t, 3)]
            @test depth(t, node_by_source_id(t, 3)) == 2
        end

        @testset "an end-to-end run on a converted tree" begin
            root = NonMarkovEvolution.BinaryNode(
                NonMarkovEvolution.NonMarkovCell(1, 0.0, 0, 0, 1.0))
            NonMarkovEvolution.left_child!(root,
                NonMarkovEvolution.NonMarkovCell(2, 1.0, 6, 6, 1.0))
            NonMarkovEvolution.right_child!(root,
                NonMarkovEvolution.NonMarkovCell(3, 1.0, 6, 6, 1.0))
            t = PhyloTree(root)
            a = toy_assembly(nchrom = 2, len = 1000)
            res = simulate_cnas(t, a, CNAModel(rate = FromEdgeMutations(0.5)); seed = 61)
            for i in 1:nnodes(t)
                @test check_invariants(profile(res, i))
            end
            @test replay(res) == [profile(res, i) for i in 1:nnodes(t)]
        end

        # A complete binary lineage tree of `depth` divisions as NonMarkovEvolution
        # BinaryNodes; `total_drivers` is kept consistent with the path from the root.
        function nme_tree(depth)
            M = NonMarkovEvolution
            root = M.BinaryNode(M.NonMarkovCell(1, 0.0, 0, 0, 1.0))
            next = Ref(2)
            frontier = [root]
            for d in 1:depth
                newf = eltype(frontier)[]
                for nd in frontier
                    tot = nd.data.total_drivers
                    M.left_child!(nd, M.NonMarkovCell(next[], Float64(d), d % 3, tot + d % 3, 1.0))
                    next[] += 1
                    M.right_child!(nd, M.NonMarkovCell(next[], d + 0.5, (d + 1) % 3, tot + (d + 1) % 3, 1.0))
                    next[] += 1
                    push!(newf, nd.left, nd.right)
                end
                frontier = newf
            end
            return root
        end

        @testset "sampling with the real sampler commutes exactly (per_node)" begin
            root = nme_tree(5)                                   # 32 leaves
            s = NonMarkovEvolution.sample_leaves(root, 6; seed = 9)
            full, sub = PhyloTree(root), PhyloTree(s)
            @test nnodes(sub) < nnodes(full)
            m = CNAModel(rate = FromEdgeMutations(1.0), target = CNWeighted(),
                         extent = ExtentMixture(p_arm = 0.2), initial = TruncalCNAs(2),
                         wgd = RateWGD(PerDivision(0.2)))
            a = toy_assembly(nchrom = 3, len = 2000)
            rf = simulate_cnas(full, a, m; seed = 62, rng_mode = :per_node)
            rs = simulate_cnas(sub, a, m; seed = 62, rng_mode = :per_node)
            for l in leaves(sub)
                sid = node(sub, l).source_id
                @test profile(rs, l) == profile(rf, node_by_source_id(full, sid))
            end
        end
    end
end
