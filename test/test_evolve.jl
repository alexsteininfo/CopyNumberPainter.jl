@testset "evolve" begin
    # Balanced 4-leaf tree, one division per edge:
    #        1
    #      /   \
    #     2     3
    #    / \   / \
    #   4   5 6   7
    bal() = phylotree([nothing, 1, 1, 2, 2, 3, 3];
                      birthtimes = [0.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0],
                      edge_divisions = [nothing, 1, 1, 1, 1, 1, 1],
                      edge_mutations = [nothing, 2, 2, 1, 1, 1, 1],
                      source_ids = collect(101:107))

    A() = toy_assembly(nchrom = 2, len = 1000)

    @testset "CNAModel defaults" begin
        m = CNAModel()
        @test m.rate isa PerDivision
        @test m.target isa CNWeighted
        @test m.extent isa ExtentMixture
        @test m.kind isa GainLoss
        @test m.wgd isa NoWGD
        @test m.viability isa RejectAndRedraw
        @test m.initial isa Diploid
    end

    @testset "zero alterations: every tip is exactly the root state" begin
        res = simulate_cnas(bal(), A(), CNAModel(rate = PerDivision(0.0)); seed = 1)
        d = diploid(A())
        for i in 1:nnodes(res.tree)
            @test profile(res, i) == d
            @test check_invariants(profile(res, i))
        end
        @test nevents(res) == 0
        @test rejection_count(res) == 0
    end

    @testset "profiles are retained for every node by default" begin
        res = simulate_cnas(bal(), A(), CNAModel(rate = PerDivision(1.0)); seed = 2)
        @test res.retain_internal
        @test all(res.profiles[i] !== nothing for i in 1:nnodes(res.tree))
        @test length(leaf_profiles(res)) == 4
        for i in 1:nnodes(res.tree)
            @test check_invariants(profile(res, i))
        end
    end

    @testset "retain_internal = false keeps only the leaves" begin
        res = simulate_cnas(bal(), A(), CNAModel(rate = PerDivision(1.0));
                            seed = 3, retain_internal = false)
        for i in internal_nodes(res.tree)
            @test res.profiles[i] === nothing
            @test_throws ArgumentError profile(res, i)
        end
        @test all(res.profiles[i] !== nothing for i in leaves(res.tree))
        @test length(leaf_profiles(res)) == 4
    end

    @testset "inheritance: a scheduled doubling reaches exactly its descendants" begin
        t = bal()
        res = simulate_cnas(t, A(),
                            CNAModel(rate = PerDivision(0.0),
                                     wgd = ScheduledWGD(2 => 1));
                            seed = 4)
        d = diploid(A())
        # node 2 and everything below it is doubled
        for i in (2, 4, 5)
            @test cn_at(slot_segments(profile(res, i), 1, 1), 1) == 2
        end
        # node 1, 3, 6, 7 are untouched
        for i in (1, 3, 6, 7)
            @test profile(res, i) == d
        end
        @test nevents(res) == 1
        @test only(res.events).node == 2
        @test only(res.events).event isa WholeGenomeDoubling
    end

    @testset "inheritance: a segmental alteration reaches exactly its descendants" begin
        t = bal()
        fixed_target = (p, rng) -> (1, 1)
        fixed_extent = (p, c, h, rng) -> (101, 200, :focal)
        # one alteration, only on the edge into node 3
        onlyon3 = CustomRate((tr, i, rng) -> i == 3 ? 1 : 0)
        res = simulate_cnas(t, A(),
                            CNAModel(rate = onlyon3, target = fixed_target,
                                     extent = fixed_extent, kind = GainLoss(1.0));
                            seed = 5)
        for i in (3, 6, 7)
            @test cn_at(slot_segments(profile(res, i), 1, 1), 150) == 2
        end
        for i in (1, 2, 4, 5)
            @test cn_at(slot_segments(profile(res, i), 1, 1), 150) == 1
        end
        @test nevents(res) == 1
        @test events_on(res, 3) |> length == 1
        @test events_on(res, 6) |> isempty
        @test length(events_below(res, 1)) == 1
        @test isempty(events_below(res, 2))
    end

    @testset "determinism" begin
        t = bal()
        m = CNAModel(rate = PerDivision(2.0), extent = ExtentMixture(p_chromosome = 0.2, p_arm = 0.2))
        a = simulate_cnas(t, A(), m; seed = 99)
        b = simulate_cnas(t, A(), m; seed = 99)
        c = simulate_cnas(t, A(), m; seed = 100)
        @test [profile(a, i) for i in 1:nnodes(t)] == [profile(b, i) for i in 1:nnodes(t)]
        @test a.events == b.events
        @test [profile(a, i) for i in 1:nnodes(t)] != [profile(c, i) for i in 1:nnodes(t)]
        @test a.seed == 99
    end

    @testset "an explicit rng is honoured and records no seed" begin
        t = bal()
        m = CNAModel(rate = PerDivision(1.0))
        a = simulate_cnas(t, A(), m; rng = Random.Xoshiro(5))
        b = simulate_cnas(t, A(), m; rng = Random.Xoshiro(5))
        @test a.events == b.events
        @test a.seed === nothing
    end

    @testset "events are stored in preorder, then by order within a node" begin
        t = bal()
        res = simulate_cnas(t, A(), CNAModel(rate = PerDivision(3.0)); seed = 6)
        rank = Dict(n => k for (k, n) in enumerate(preorder(t)))
        keys_ = [(rank[e.node], e.order) for e in res.events]
        @test issorted(keys_)
        for i in 1:nnodes(t)
            ons = events_on(res, i)
            @test [e.order for e in ons] == collect(1:length(ons))
        end
    end

    @testset "TruncalCNAs sets the root state and logs its events there" begin
        t = bal()
        res = simulate_cnas(t, A(),
                            CNAModel(rate = PerDivision(0.0),
                                     initial = TruncalCNAs(3),
                                     kind = GainLoss(1.0));
                            seed = 7)
        @test nevents(res) == 3
        @test all(e.node == treeroot(t) for e in res.events)
        @test [e.order for e in res.events] == [1, 2, 3]
        @test profile(res, treeroot(t)) != diploid(A())
        # with no further alterations every tip equals the root
        for i in leaves(t)
            @test profile(res, i) == profile(res, treeroot(t))
        end
        @test_throws ArgumentError TruncalCNAs(-1)
    end

    @testset "Given sets the root state and must match the assembly" begin
        t = bal()
        start = diploid(A())
        apply!(start, SegmentalCNA(1, 1, 1, 1000, 1, :chromosome))
        res = simulate_cnas(t, A(), CNAModel(rate = PerDivision(0.0), initial = Given(start)); seed = 8)
        @test profile(res, treeroot(t)) == start
        for i in leaves(t)
            @test profile(res, i) == start
        end
        # the given profile is not mutated by the simulation
        @test cn_at(slot_segments(start, 1, 1), 1) == 2
        @test_throws ArgumentError simulate_cnas(t, toy_assembly(nchrom = 3, len = 1000),
                                                 CNAModel(initial = Given(start)); seed = 9)
    end

    @testset "event-log replay reproduces every profile" begin
        t = bal()
        m = CNAModel(rate = PerDivision(2.0), initial = TruncalCNAs(2),
                     wgd = ScheduledWGD(3 => 1),
                     extent = ExtentMixture(p_chromosome = 0.1, p_arm = 0.2))
        res = simulate_cnas(t, A(), m; seed = 10)
        rp = replay(res)
        for i in 1:nnodes(t)
            @test rp[i] == profile(res, i)
        end
        # replay also works when nothing was retained
        lean = simulate_cnas(t, A(), m; seed = 10, retain_internal = false)
        rp2 = replay(lean)
        for i in 1:nnodes(t)
            @test rp2[i] == rp[i]
        end
    end

    @testset "the result carries its own root state and event index" begin
        t = bal()
        m = CNAModel(rate = PerDivision(3.0), initial = TruncalCNAs(3; wgd = 1))
        res = simulate_cnas(t, A(), m; seed = 19)
        @test res.root_base == diploid(A())                 # before the truncal events
        @test sum(length, res.event_ranges) == nevents(res)
        for i in 1:nnodes(t)
            @test all(e.node == i for e in events_on(res, i))
            @test events_on(res, i) == [e for e in res.events if e.node == i]
        end
        below(i) = Set(j for j in preorder(t) if j != i && i in ancestors(t, j))
        for i in 1:nnodes(t)
            @test events_below(res, i) == [e for e in res.events if e.node in below(i)]
        end
        # replay no longer consults the model
        stripped = CopyNumberEvolution.CNAEvolution(res.tree, res.assembly, nothing,
            res.profiles, res.events, res.event_ranges, res.root_base, res.rejections,
            res.seed, res.rng_mode, res.retain_internal)
        @test replay(stripped) == [profile(res, i) for i in 1:nnodes(t)]
    end

    @testset "_event_ranges rejects logs that replay could not trust" begin
        ev(node, order) = LoggedEvent(node, order, WholeGenomeDoubling(:multiply))
        er = CopyNumberEvolution._event_ranges
        @test er([ev(1, 1), ev(2, 1), ev(2, 2), ev(3, 1)], 4) == [1:1, 2:3, 4:4, 1:0]
        @test_throws ArgumentError er([ev(2, 1), ev(3, 1), ev(2, 2)], 3)   # interleaved
        @test_throws ArgumentError er([ev(1, 1), ev(5, 1)], 3)             # node out of range
        @test_throws ArgumentError er([ev(2, 2), ev(2, 1)], 3)             # orders [2, 1]
    end

    @testset "proposals that would change nothing are redrawn, not logged" begin
        a = toy_assembly(nchrom = 1, len = 1000)
        start = diploid(a)
        apply!(start, SegmentalCNA(1, 1, 1, 600, -1, :focal))     # hap1 1:600 absent
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        m = CNAModel(rate = CustomRate((tr, i, rng) -> 50), target = (p, rng) -> (1, 1),
                     extent = ExtentMixture(lengthdist = Distributions.Uniform(10.0, 100.0)),
                     viability = AllowAll(), initial = Given(start))
        res = simulate_cnas(t, a, m; seed = 48)
        @test nevents(res) == 50
        p = copy(start)
        for le in events_on(res, 2)
            before = copy(p)
            apply!(p, le.event)
            @test p != before                           # every logged event did something
        end
        @test get(res.rejections, :no_effect, 0) > 0
    end

    @testset "a genome with no material left is an error, not a hang" begin
        a = toy_assembly(nchrom = 1, len = 100)
        dead = diploid(a)
        apply!(dead, SegmentalCNA(1, 1, 1, 100, -1, :chromosome))
        apply!(dead, SegmentalCNA(1, 2, 1, 100, -1, :chromosome))
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        m = CNAModel(rate = CustomRate((tr, i, rng) -> 1), target = UniformChromosome(),
                     viability = AllowAll(), initial = Given(dead))
        err = try simulate_cnas(t, a, m; seed = 1); nothing catch e; e end
        @test err isa ErrorException && occursin("no copy-number material", err.msg)
    end

    @testset "rejections are tallied by reason" begin
        t = bal()
        # mostly losses, on a male karyotype: hemizygous chrX and chrY losses get rejected
        # UniformChromosome is pinned: CNWeighted never targets deleted material.
        # Gains are mixed in because losses on already-deleted slots are redrawn as
        # :no_effect, so pure losses would eventually leave only rejectable ones and
        # exhaust max_attempts instead of tallying.
        m = CNAModel(rate = PerDivision(4.0), target = UniformChromosome(),
                     kind = GainLoss(0.3),
                     extent = ExtentMixture(p_chromosome = 1.0),
                     viability = RejectAndRedraw(min_total_cn = 1, max_attempts = 10_000))
        res = simulate_cnas(t, toy_sex_assembly(:male), m; seed = 11)
        @test rejection_count(res) > 0
        @test haskey(res.rejections, :min_total_cn)
        @test res.rejections[:min_total_cn] + get(res.rejections, :no_effect, 0) ==
              rejection_count(res)
    end

    @testset "exhausting max_attempts throws with actionable advice" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        # a single hemizygous chromosome, whole-chromosome losses only: every
        # proposal drives total copy number to 0 and is rejected
        m = CNAModel(rate = PerDivision(5.0), kind = GainLoss(0.0),
                     extent = ExtentMixture(p_chromosome = 1.0),
                     viability = RejectAndRedraw(min_total_cn = 1, max_attempts = 5))
        err = try
            simulate_cnas(t, hemizygous_assembly(), m; seed = 12)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("max_attempts", err.msg)
        @test occursin("min_total_cn", err.msg)
    end

    @testset "rng_mode = :per_node stays reproducible" begin
        t = bal()
        m = CNAModel(rate = PerDivision(2.0))
        @test_throws ArgumentError simulate_cnas(t, A(), m; seed = 1, rng_mode = :bogus)
        @test_throws ArgumentError simulate_cnas(t, A(), m; seed = UInt64(typemax(Int)) + 1)
        a = simulate_cnas(t, A(), m; seed = 13, rng_mode = :per_node)
        b = simulate_cnas(t, A(), m; seed = 13, rng_mode = :per_node)
        @test a.events == b.events
        for i in 1:nnodes(t)
            @test check_invariants(profile(a, i))
        end
    end

    @testset "seeded streams are pinned across Julia versions" begin
        r = CopyNumberEvolution._stable_rng(UInt64(1))
        @test rand(r, UInt64) == 0x6bfbe5ada17babc0
        @test rand(r, UInt64) == 0xd9791c54a38dd5f4
        @test rand(CopyNumberEvolution._stable_rng(UInt64(1), UInt64(2), UInt64(7)), UInt64) ==
              0x336c294898a4faba
    end

    @testset "a run without seed or rng draws and records a seed" begin
        Random.seed!(123); a = simulate_cnas(bal(), A(), CNAModel(rate = PerDivision(2.0)))
        Random.seed!(123); b = simulate_cnas(bal(), A(), CNAModel(rate = PerDivision(2.0)))
        @test a.seed isa Int && a.seed == b.seed && a.events == b.events
        c = simulate_cnas(bal(), A(), CNAModel(rate = PerDivision(2.0)); seed = a.seed)
        @test c.events == a.events                        # reproducible from the record alone
        @test a.rng_mode === :global
        @test simulate_cnas(bal(), A(), CNAModel(); rng = Random.Xoshiro(1)).seed === nothing
        @test_throws ArgumentError simulate_cnas(bal(), A(), CNAModel();
                                                 rng = Random.Xoshiro(1), rng_mode = :per_node)
        @test simulate_cnas(bal(), A(), CNAModel(); rng_mode = :per_node).rng_mode === :per_node
    end

    @testset "per_node keys: a missing source_id never collides with a present one" begin
        # node 2 has no source_id (falls back to dense id 2); node 3's source_id is 2
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                      source_ids = [nothing, nothing, 2])
        res = simulate_cnas(t, A(), CNAModel(rate = PerDivision(5.0)); seed = 3, rng_mode = :per_node)
        @test [e.event for e in events_on(res, 2)] != [e.event for e in events_on(res, 3)]
    end

    @testset "a unary chain draws once per edge" begin
        t = phylotree([nothing, 1, 2, 3]; edge_divisions = [nothing, 1, 1, 1])
        onlyone = CNAModel(rate = PerDivision(0.0), wgd = ScheduledWGD(2 => 1, 3 => 1, 4 => 1))
        res = simulate_cnas(t, A(), onlyone; seed = 14)
        @test nevents(res) == 3
        @test cn_at(slot_segments(profile(res, 4), 1, 1), 1) == 8   # doubled three times
        @test cn_at(slot_segments(profile(res, 1), 1, 1), 1) == 1
    end

    @testset "deep trees do not overflow the stack" begin
        n = 20_000
        parents = Vector{Union{Int,Nothing}}(undef, n)
        parents[1] = nothing
        for i in 2:n
            parents[i] = i - 1
        end
        t = phylotree(parents; edge_divisions = vcat(nothing, fill(1, n - 1)))
        res = simulate_cnas(t, A(), CNAModel(rate = PerDivision(0.05));
                            seed = 15, retain_internal = false)
        @test check_invariants(profile(res, n))
        @test nevents(res) > 0
    end

    @testset "show methods do not error" begin
        res = simulate_cnas(bal(), A(), CNAModel(rate = PerDivision(1.0)); seed = 16)
        @test occursin("CNAEvolution", sprint(show, res))
        @test occursin("CNAModel", sprint(show, res.model))
    end


    @testset "Given rejects a non-canonical profile" begin
        p = diploid(A())
        p.segments[1] = [CopyNumberEvolution.Segment(1, 500, 1), CopyNumberEvolution.Segment(501, 1000, 1)]
        @test_throws ErrorException Given(p)
    end

    @testset "CNAModel checks component types at construction" begin
        @test_throws ArgumentError CNAModel(rate = 0.5)
        @test_throws ArgumentError CNAModel(wgd = 1)
        @test_throws ArgumentError CNAModel(viability = :none)
        @test_throws ArgumentError CNAModel(initial = diploid(A()))
        # the three draws still accept plain functions
        @test CNAModel(target = (p, rng) -> (1, 1)) isa CNAModel
    end

    @testset "a custom InitialState needs only initial_profile" begin
        t = bal()
        res = simulate_cnas(t, A(), CNAModel(rate = PerDivision(1.0), initial = DoubledChr1()); seed = 18)
        @test cn_at(slot_segments(profile(res, treeroot(t)), 1, 1), 1) == 2
        @test replay(res) == [profile(res, i) for i in 1:nnodes(t)]
        @test initial_profile(TruncalCNAs(3), A()) == diploid(A())
        @test initial_profile(Diploid(), A()) == diploid(A())
    end

    @testset "a doubling lands at a uniformly random position on its edge" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        m = CNAModel(rate = CustomRate((tr, i, rng) -> 9), wgd = ScheduledWGD(2 => 1),
                     viability = AllowAll())
        pos = Int[]
        for s in 1:600
            res = simulate_cnas(t, A(), m; seed = s)
            evs = events_on(res, 2)
            @test length(evs) == 10
            @test [e.order for e in evs] == 1:10
            push!(pos, only(e.order for e in evs if e.event isa WholeGenomeDoubling))
            s <= 20 && @test replay(res)[2] == profile(res, 2)
        end
        @test Set(pos) == Set(1:10)
        @test sum(pos) / length(pos) ≈ 5.5 atol = 0.4
    end

    @testset "adding a doubling leaves other edges' draws unchanged (per_node)" begin
        t = bal()
        a = simulate_cnas(t, A(), CNAModel(rate = PerDivision(2.0)); seed = 17, rng_mode = :per_node)
        b = simulate_cnas(t, A(), CNAModel(rate = PerDivision(2.0), wgd = ScheduledWGD(3 => 1));
                          seed = 17, rng_mode = :per_node)
        for i in (1, 2, 4, 5)                    # outside the doubled subtree
            @test [e.event for e in events_on(a, i)] == [e.event for e in events_on(b, i)]
        end
    end

    @testset "TruncalCNAs can place doublings among the truncal events" begin
        t = bal()
        m = CNAModel(rate = PerDivision(0.0), initial = TruncalCNAs(4; wgd = 1),
                     wgd = ScheduledWGD(2 => 1; mode = :increment))
        seen = Set{Int}()
        for s in 1:200
            res = simulate_cnas(t, A(), m; seed = s)
            root = events_on(res, treeroot(t))
            @test length(root) == 5
            d = only(e for e in root if e.event isa WholeGenomeDoubling)
            @test d.event.mode === :increment               # inherited from the policy
            push!(seen, d.order)
            s <= 10 && @test replay(res) == [profile(res, i) for i in 1:nnodes(t)]
        end
        @test seen == Set(1:5)
        # doublings only, and more than one
        r0 = simulate_cnas(t, A(), CNAModel(rate = PerDivision(0.0), initial = TruncalCNAs(0; wgd = 2)); seed = 1)
        @test nevents(r0) == 2
        @test cn_at(slot_segments(profile(r0, treeroot(t)), 1, 1), 1) == 4
        @test replay(r0) == [profile(r0, i) for i in 1:nnodes(t)]
        @test TruncalCNAs(2; wgd = 1, mode = :multiply).mode === :multiply
        @test TruncalCNAs(3).wgd == 0
        @test_throws ArgumentError TruncalCNAs(2; wgd = -1)
        @test_throws ArgumentError TruncalCNAs(2; mode = :bogus)
    end
    @testset "ExactlyNWGD under per_node warns that it does not commute" begin
        @test_logs (:warn, r"ExactlyNWGD") simulate_cnas(bal(), A(),
            CNAModel(wgd = ExactlyNWGD(1)); seed = 1, rng_mode = :per_node)
        @test_logs simulate_cnas(bal(), A(), CNAModel(wgd = ExactlyNWGD(1)); seed = 1)
    end

end
