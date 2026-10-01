using TOML

@testset "records" begin
    @testset "model descriptions are TOML-safe, closures included" begin
        m = CNAModel(rate = CustomRate((t, i, rng) -> 1), target = (p, rng) -> (1, 1),
                     extent = ExtentMixture(p_arm = 0.2),
                     wgd = ScheduledWGD(2 => 1; mode = :increment),
                     viability = AllRules([RejectAndRedraw(), AllowAll()]),
                     initial = TruncalCNAs(3; wgd = 1))
        d = CopyNumberEvolution._describe(m)
        buf = IOBuffer()
        TOML.print(buf, Dict("model" => d))                 # must not throw
        back = TOML.parse(String(take!(buf)))["model"]
        @test back["type"] == "CNAModel"
        @test back["rate"]["type"] == "CustomRate"
        @test back["rate"]["f"]["type"] == "function"
        @test back["target"]["type"] == "function"
        @test back["extent"]["p_arm"] == 0.2  # TOML preserves numeric type
        @test back["extent"]["lengthdist"]["type"] == "LogUniform"
        @test back["wgd"]["at"] == Dict("2" => 1) && back["wgd"]["mode"] == "increment"
        @test length(back["viability"]["rules"]) == 2
        @test back["initial"]["mode"] == "nothing"
        @test back["initial"]["type"] == "TruncalCNAs"
    end

    @testset "assemblies round-trip through TOML" begin
        for a in (hg38(:male), hg19(:female), toy_sex_assembly(:male), hemizygous_assembly())
            buf = IOBuffer()
            TOML.print(buf, CopyNumberEvolution._assembly_meta(a))
            b = CopyNumberEvolution._assembly_from_meta(TOML.parse(String(take!(buf))))
            @test same_assembly(a, b)
        end
    end

    @testset "ModelRecord shows its original repr" begin
        r = ModelRecord(Dict{String,Any}("type" => "CNAModel"), "CNAModel(rate=…)")
        @test sprint(show, r) == "CNAModel(rate=…)"
        @test CopyNumberEvolution._describe(r) == r.description
    end

    @testset "_describe handles cyclic user structs without stack overflow" begin
        cyc = CyclicUserStruct(nothing)
        cyc.self = cyc
        d = CopyNumberEvolution._describe(cyc)
        buf = IOBuffer()
        # Should not throw StackOverflowError; depth limit prevents infinite recursion
        TOML.print(buf, Dict("cyc" => d))
        back = TOML.parse(String(take!(buf)))["cyc"]
        @test back["type"] == "CyclicUserStruct"
    end

    @testset "_describe handles structs with type-valued fields" begin
        ts = TypeFieldStruct(Int, 42)
        d = CopyNumberEvolution._describe(ts)
        buf = IOBuffer()
        TOML.print(buf, Dict("ts" => d))
        back = TOML.parse(String(take!(buf)))["ts"]
        @test back["type_field"] == "Int64"
        @test back["value"] == 42  # TOML preserves numeric type
    end

    @testset "_describe handles matrix fields with summary" begin
        mfs = MatrixFieldStruct(rand(Int, 3, 4))
        d = CopyNumberEvolution._describe(mfs)
        buf = IOBuffer()
        TOML.print(buf, Dict("mfs" => d))
        back = TOML.parse(String(take!(buf)))["mfs"]
        @test back["type"] == "MatrixFieldStruct"
        # Matrix should appear as its summary
        @test contains(back["mat"], "Matrix")
    end

    @testset "_describe handles long vector fields with summary" begin
        lvfs = LongVectorFieldStruct(collect(1:2000))
        d = CopyNumberEvolution._describe(lvfs)
        buf = IOBuffer()
        TOML.print(buf, Dict("lvfs" => d))
        back = TOML.parse(String(take!(buf)))["lvfs"]
        @test back["type"] == "LongVectorFieldStruct"
        # Long vector should appear as its summary
        @test contains(back["vec"], "Vector") || contains(back["vec"], "Array")
    end

    @testset "_describe handles Tuple fields as vectors" begin
        tfs = TupleFieldStruct((1, "hello", 3.14))
        d = CopyNumberEvolution._describe(tfs)
        buf = IOBuffer()
        TOML.print(buf, Dict("tfs" => d))
        back = TOML.parse(String(take!(buf)))["tfs"]
        @test back["type"] == "TupleFieldStruct"
        # Tuple should appear as a vector of described elements
        @test isa(back["tup"], Vector)
        @test back["tup"][1] == 1
        @test back["tup"][2] == "hello"
        @test back["tup"][3] == 3.14
    end

    @testset "_describe handles Set fields as sorted vectors" begin
        sfs = SetFieldStruct(Set([3, 1, 2]))
        d = CopyNumberEvolution._describe(sfs)
        buf = IOBuffer()
        TOML.print(buf, Dict("sfs" => d))
        back = TOML.parse(String(take!(buf)))["sfs"]
        @test back["type"] == "SetFieldStruct"
        # Set should appear as a sorted vector of described elements
        @test isa(back["s"], Vector)
        @test back["s"] == [1, 2, 3]
    end

    @testset "_describe handles self-containing vectors without stack overflow" begin
        scvs = SelfContainingVectorStruct(Any[])
        push!(scvs.vec, scvs)
        d = CopyNumberEvolution._describe(scvs)
        buf = IOBuffer()
        # Should not throw StackOverflowError; depth limit in container methods prevents it
        TOML.print(buf, Dict("scvs" => d))
        back = TOML.parse(String(take!(buf)))["scvs"]
        @test back["type"] == "SelfContainingVectorStruct"
        # The vector field should be present (depth limit prevents infinite recursion)
        @test haskey(back, "vec")
    end

    function sim()
        t = binary_lineage(3)
        m = CNAModel(rate = PerDivision(1.5), target = CNWeighted(),
                     extent = ExtentMixture(p_arm = 0.2), initial = TruncalCNAs(2; wgd = 1))
        return simulate_cnas(t, toy_sex_assembly(:male, len = 2000), m; seed = 52)
    end
    allprofiles(r) = [profile(r, i) for i in 1:nnodes(r.tree)]

    @testset "a simulation bundle round-trips" begin
        res = sim()
        for compress in (false, true), profiles in (:all, :leaves, :none)
            prefix = joinpath(mktempdir(), "run")
            save_simulation(prefix, res; profiles = profiles, compress = compress)
            back = load_simulation(prefix)
            @test back.model isa ModelRecord
            @test back.events == res.events
            @test allprofiles(back) == allprofiles(res)
            @test back.root_base == res.root_base
            @test back.seed == res.seed && back.rng_mode === res.rng_mode
            @test back.rejections == res.rejections
            @test same_assembly(back.assembly, res.assembly)
            @test replay(back) == allprofiles(res)
            meta = TOML.parsefile(prefix * "_meta.toml")
            @test meta["format_version"] == 1
            @test meta["package_version"] == string(pkgversion(CopyNumberEvolution))
            @test meta["newick_branchlength"] == "divisions"
            @test isfile(prefix * ".nwk")
            save_simulation(prefix * "2", back)             # a loaded result saves again
            @test load_simulation(prefix * "2").events == res.events
        end
    end

    @testset "a lean run saved with leaves only loads complete" begin
        t = binary_lineage(3)
        m = CNAModel(rate = PerDivision(2.0))
        a = toy_assembly(nchrom = 2, len = 1000)
        full = simulate_cnas(t, a, m; seed = 53)
        lean = simulate_cnas(t, a, m; seed = 53, retain_internal = false)
        prefix = joinpath(mktempdir(), "lean")
        save_simulation(prefix, lean; profiles = :leaves)
        @test allprofiles(load_simulation(prefix)) == allprofiles(full)
    end

    @testset "a missing optional part warns and is rebuilt" begin
        res = sim()
        prefix = joinpath(mktempdir(), "run")
        save_simulation(prefix, res)
        rm(prefix * "_profiles.tsv")
        back = @test_logs (:warn, r"missing") load_simulation(prefix)
        @test allprofiles(back) == allprofiles(res)
        rm(prefix * "_events.tsv")
        @test_throws Exception load_simulation(prefix)       # the log is required
    end

    @testset "a newer bundle format is refused" begin
        res = sim()
        prefix = joinpath(mktempdir(), "run")
        save_simulation(prefix, res)
        meta = TOML.parsefile(prefix * "_meta.toml")
        meta["format_version"] = 99
        open(io -> TOML.print(io, meta), prefix * "_meta.toml", "w")
        @test_throws ArgumentError load_simulation(prefix)
    end

    @testset "a matrix bundle round-trips, mask and phasing included" begin
        res = sim()
        a = res.assembly
        g = BinGrid(a, 300; mask = vcat(centromere_mask(a), ["chr1" => 1:450]))
        for m in (CNMatrix(res, g; rule = AreaWeightedMean()),
                  major_minor(CNMatrix(res, g)),
                  CNMatrix(res, g; allele = false)), compress in (false, true)
            prefix = joinpath(mktempdir(), "mat")
            save_matrix(prefix, m; compress = compress)
            back = load_matrix(prefix)
            @test back.grid.bins == g.bins && back.grid.mask == g.mask
            @test back.cells == m.cells && back.names == m.names
            @test back.total == m.total && back.allele == m.allele
            @test back.phasing === m.phasing
            @test typeof(back.rule) == typeof(m.rule)
        end
    end

    @testset "seeds of any size survive a save/load round trip" begin
        t = binary_lineage(3)
        a = toy_assembly(nchrom = 2, len = 1000)
        m = CNAModel(rate = PerDivision(1.0))
        # typemax(Int) is an Int128 once TOML has parsed it back
        big = simulate_cnas(t, a, m; seed = typemax(Int))
        prefix = joinpath(mktempdir(), "big")
        save_simulation(prefix, big)
        back = load_simulation(prefix)
        @test back.seed === typemax(Int) && back.seed isa Int
        @test back.events == big.events
        # a drawn seed is uniform on 0:typemax(Int), so half of them exceed 2^62
        for s in (2^62, 2^62 + 1)
            r = simulate_cnas(t, a, m; seed = s)
            p = joinpath(mktempdir(), "s")
            save_simulation(p, r)
            @test load_simulation(p).seed === s
        end
        drawn = simulate_cnas(t, a, m)
        pd = joinpath(mktempdir(), "drawn")
        save_simulation(pd, drawn)
        @test load_simulation(pd).seed === drawn.seed
    end

    @testset "a Given initial state is described briefly" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0), initial = Given(diploid(a)));
                            seed = 3)
        prefix = joinpath(mktempdir(), "given")
        save_simulation(prefix, res)
        meta = TOML.parsefile(prefix * "_meta.toml")
        @test meta["model"]["initial"]["profile"] == "CNProfile(toy, 2 segments)"
        @test load_simulation(prefix).events == res.events
    end

    @testset "a refused name leaves no partial bundle" begin
        dup = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                        labels = ["R", "x", "x"])
        res = simulate_cnas(dup, toy_assembly(nchrom = 1, len = 400),
                            CNAModel(rate = PerDivision(1.0)); seed = 4)
        dir = mktempdir()
        @test_throws ArgumentError save_simulation(joinpath(dir, "run"), res)
        @test isempty(readdir(dir))
        na = phylotree([nothing, 1]; edge_divisions = [nothing, 1], labels = ["R", "NA"])
        resna = simulate_cnas(na, toy_assembly(nchrom = 1, len = 400),
                              CNAModel(rate = PerDivision(1.0)); seed = 4)
        @test_throws ArgumentError save_simulation(joinpath(dir, "run"), resna)
        @test isempty(readdir(dir))
    end
end
