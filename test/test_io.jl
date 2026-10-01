@testset "io" begin
    setup() = begin
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                      labels = [nothing, "left", "right"])
        a = toy_sex_assembly(:male, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(1.0),
                                           extent = ExtentMixture(p_chromosome = 0.5));
                            seed = 41)
        (t, a, res, BinGrid(a, 100))
    end

    readtsv(path) = [split(l, '\t') for l in readlines(path)]

    @testset "write_profiles" begin
        t, a, res, g = setup()
        path = joinpath(mktempdir(), "profiles.tsv")
        write_profiles(path, res)
        rows = readtsv(path)
        @test rows[1] == ["node_id", "name", "chrom", "haplotype", "start", "stop", "cn"]
        @test length(rows) > 1
        body = rows[2:end]
        @test all(length(r) == 7 for r in body)
        # every chromosome name present is a real one, and coordinates are 1-based
        @test all(r[3] in [chromname(a, c) for c in 1:nchromosomes(a)] for r in body)
        @test minimum(parse(Int, r[5]) for r in body) == 1
        # every chromosome written has slots — the set of names is exactly the
        # positive-ploidy ones, no more and no less
        written = Set(r[3] for r in body)
        @test written == Set(chromname(a, c) for c in 1:nchromosomes(a) if ploidy(a, c) > 0)
    end

    @testset "a zero-ploidy chromosome never appears in a profile table" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_sex_assembly(:female, len = 400)      # chrY has ploidy 0 here
        @test ploidy(a, chromindex(a, "chrY")) == 0
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0)); seed = 45)
        path = joinpath(mktempdir(), "female.tsv")
        write_profiles(path, res)
        chroms = Set(r[3] for r in readtsv(path)[2:end])
        @test !("chrY" in chroms)
        @test "chrX" in chroms
        @test "chr1" in chroms
    end

    @testset "write_profiles honours the cells argument" begin
        t, a, res, g = setup()
        path = joinpath(mktempdir(), "leaves.tsv")
        write_profiles(path, res; cells = leaves(t))
        ids = Set(parse(Int, r[1]) for r in readtsv(path)[2:end])
        @test ids == Set(leaves(t))
    end

    @testset "write_events" begin
        t, a, res, g = setup()
        path = joinpath(mktempdir(), "events.tsv")
        write_events(path, res)
        rows = readtsv(path)
        @test rows[1] == ["node_id", "name", "order", "type", "chrom", "haplotype",
                          "start", "stop", "delta", "scale", "mode"]
        @test length(rows) - 1 == nevents(res)
        for r in rows[2:end]
            @test r[4] in ("segmental", "wgd")
            if r[4] == "segmental"
                @test r[11] == "NA"
                @test parse(Int, r[9]) != 0
            else
                @test r[5] == "NA" && r[11] in ("multiply", "increment")
            end
        end
    end

    @testset "write_events on a run with a doubling records the mode" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0),
                                           wgd = ScheduledWGD(2 => 1; mode = :increment));
                            seed = 42)
        path = joinpath(mktempdir(), "e.tsv")
        write_events(path, res)
        rows = readtsv(path)
        @test length(rows) == 2
        @test rows[2][4] == "wgd"
        @test rows[2][11] == "increment"
    end

    @testset "write_bins" begin
        t, a, res, g = setup()
        path = joinpath(mktempdir(), "bins.tsv")
        write_bins(path, g)
        rows = readtsv(path)
        @test rows[1] == ["bin_index", "chrom", "start", "stop"]
        @test length(rows) - 1 == nbins(g)
        @test rows[2] == ["1", "chr1", "1", "100"]
    end

    @testset "write_medicc2 shape and coordinates" begin
        t, a, res, g = setup()
        m = CNMatrix(res, g)
        path = joinpath(mktempdir(), "medicc.tsv")
        write_medicc2(path, m)
        rows = readtsv(path)
        @test rows[1] == ["sample_id", "chrom", "start", "end", "cn_a", "cn_b"]
        autosomal = sum(length(bins_of(g, c)) for c in autosomes(a))
        @test length(rows) - 1 == (ncells(m) + 1) * autosomal
        # BED convention: 0-based start, exclusive end
        @test rows[2][2:4] == ["chr1", "0", "100"]
        # the reference rows come first and are all 1/1
        diploid_rows = [r for r in rows[2:end] if r[1] == "diploid"]
        @test length(diploid_rows) == autosomal
        @test all(r[5] == "1" && r[6] == "1" for r in diploid_rows)
        # sample ids match cellname, so they agree with an exported newick
        @test Set(r[1] for r in rows[2:end]) == Set(vcat("diploid", m.names))
    end

    @testset "write_medicc2 excludes sex chromosomes by default" begin
        t, a, res, g = setup()
        m = CNMatrix(res, g)
        path = joinpath(mktempdir(), "auto.tsv")
        write_medicc2(path, m)
        chroms = Set(r[2] for r in readtsv(path)[2:end])
        @test !("chrX" in chroms) && !("chrY" in chroms)
        path2 = joinpath(mktempdir(), "withxy.tsv")
        write_medicc2(path2, m; include_xy = true)
        chroms2 = Set(r[2] for r in readtsv(path2)[2:end])
        @test "chrX" in chroms2 && "chrY" in chroms2
    end

    @testset "write_medicc2 warns above MEDICC2's copy-number ceiling" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        start = diploid(a)
        apply!(start, SegmentalCNA(1, 1, 1, 400, 11, :chromosome))   # cn 12
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0), initial = Given(start));
                            seed = 43)
        m = CNMatrix(res, BinGrid(a, 100))
        @test max_cn(m) > 8
        path = joinpath(mktempdir(), "high.tsv")
        @test_logs (:warn, r"MEDICC2") write_medicc2(path, m)
    end

    @testset "write_medicc2 does not warn when only the total exceeds 8" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        start = diploid(a)
        apply!(start, SegmentalCNA(1, 1, 1, 400, 5, :chromosome))   # major 6
        apply!(start, SegmentalCNA(1, 2, 1, 400, 2, :chromosome))   # minor 3
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0), initial = Given(start));
                            seed = 45)
        m = CNMatrix(res, BinGrid(a, 100))
        @test max_cn(m) == 9
        @test_logs write_medicc2(joinpath(mktempdir(), "total9.tsv"), m)
    end

    @testset "write_medicc2 needs allele tracks" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0)); seed = 44)
        m = CNMatrix(res, BinGrid(a, 100); allele = false)
        @test_throws ArgumentError write_medicc2(joinpath(mktempdir(), "x.tsv"), m)
    end

    @testset "an exported newick and matrix agree on names" begin
        t, a, res, g = setup()
        m = CNMatrix(res, g)
        dir = mktempdir()
        write_newick(joinpath(dir, "t.nwk"), t; branchlength = :divisions)
        write_medicc2(joinpath(dir, "m.tsv"), m)
        nwk = read(joinpath(dir, "t.nwk"), String)
        for nm in m.names
            @test occursin(nm, nwk)
        end
    end

    @testset "a .gz path writes a complete gzip stream" begin
        t, a, res, g = setup()
        dir = mktempdir()
        for (w, x) in ((write_events, res), (write_profiles, res), (write_bins, g),
                       (write_medicc2, CNMatrix(res, g)))
            plain, gz = joinpath(dir, "x.tsv"), joinpath(dir, "x.tsv.gz")
            w(plain, x); w(gz, x)
            @test read(gz)[1:2] == [0x1f, 0x8b]                      # gzip magic bytes
            @test read(CopyNumberPainter.CodecZlib.GzipDecompressorStream(open(gz)), String) ==
                  read(plain, String)
        end
    end

    @testset "writers accept an IO as well as a path" begin
        t, a, res, g = setup()
        buf = IOBuffer()
        write_bins(buf, g)
        @test occursin("bin_index", String(take!(buf)))
    end

    @testset "the MEDICC2 normal follows the assembly's ploidy" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        for (sex, x_b, has_y) in ((:male, "0", true), (:female, "1", false))
            a = toy_sex_assembly(sex, len = 400)
            res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0)); seed = 46)
            buf = IOBuffer()
            write_medicc2(buf, CNMatrix(res, BinGrid(a, 100)); include_xy = true)
            rows = [split(l, '\t') for l in split(String(take!(buf)), '\n'; keepempty = false)][2:end]
            normal = [r for r in rows if r[1] == "diploid"]
            cell = [r for r in rows if r[1] != "diploid"]
            @test all(r[5:6] == ["1", x_b] for r in normal if r[2] == "chrX")
            @test all(r[5:6] == ["1", "1"] for r in normal if r[2] == "chr1")
            @test any(r[2] == "chrY" for r in normal) == has_y
            @test all(r[5:6] == ["1", "0"] for r in normal if r[2] == "chrY")
            # an unaltered cell is identical to the normal, bin for bin
            @test [r[2:6] for r in cell] == [r[2:6] for r in normal]
        end
    end

    @testset "write_medicc2 refuses more than two haplotype tracks" begin
        spec = [CopyNumberPainter.ChromosomeSpec("chr1", 400, 161:240)]
        a = CopyNumberPainter.GenomeAssembly("tri", :female, spec, [3])
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0)); seed = 49)
        @test_throws ArgumentError write_medicc2(IOBuffer(), CNMatrix(res, BinGrid(a, 100)))
    end

    @testset "write_profiles after a lean run writes the same file" begin
        t = phylotree([nothing, 1, 1, 2, 2]; edge_divisions = [nothing, 1, 1, 1, 1])
        a = toy_assembly(nchrom = 2, len = 1000)
        m = CNAModel(rate = PerDivision(2.0))
        full = simulate_cnas(t, a, m; seed = 47)
        lean = simulate_cnas(t, a, m; seed = 47, retain_internal = false)
        b1, b2 = IOBuffer(), IOBuffer()
        write_profiles(b1, full)
        write_profiles(b2, lean)
        @test String(take!(b1)) == String(take!(b2))
    end

    @testset "major/minor export keeps a hemizygous chrX as 1/0" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_sex_assembly(:male, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0)); seed = 50)
        buf = IOBuffer()
        write_medicc2(buf, major_minor(CNMatrix(res, BinGrid(a, 100))); include_xy = true)
        rows = [split(l, '\t') for l in split(String(take!(buf)), '\n'; keepempty = false)][2:end]
        @test all(r[5:6] == ["1", "0"] for r in rows if r[2] == "chrX")   # normal and cell
    end

    @testset "tables read back exactly" begin
        t = phylotree([nothing, 1, 1, 2, 2];
                      birthtimes = [0.0, 1.25, 2.5, 3.0, 4.75],
                      edge_divisions = [nothing, 2, 3, 1, 4],
                      edge_mutations = [nothing, 7, 0, 5, 11],
                      labels = ["R", nothing, "c", "a", "b"],
                      source_ids = [10, 20, 30, 40, 50])
        a = toy_sex_assembly(:male, len = 1000)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(3.0),
                            extent = ExtentMixture(p_arm = 0.3),
                            wgd = ScheduledWGD(3 => 1; mode = :increment)); seed = 51)
        dir = mktempdir()
        for ext in (".tsv", ".tsv.gz")
            p(x) = joinpath(dir, x * ext)
            write_tree(p("tree"), t)
            rt = read_tree(p("tree"))
            fields(tr) = [(n.parent, n.children, n.birthtime, n.edge_divisions,
                           n.edge_mutations, n.label, n.source_id) for n in tr.nodes]
            @test fields(rt) == fields(t)
            write_events(p("events"), res)
            @test read_events(p("events"), a) == res.events
            write_profiles(p("profiles"), res)
            @test read_profiles(p("profiles"), a) == Dict(i => profile(res, i) for i in 1:nnodes(t))
            g = BinGrid(a, 300; mask = centromere_mask(a))
            write_bins(p("bins"), g)
            @test read_bins(p("bins"), a) == g.bins
        end
    end

    @testset "readers refuse malformed tables" begin
        a = toy_assembly(nchrom = 1, len = 100)
        @test_throws ArgumentError read_events(IOBuffer("wrong\theader\n"), a)
        bad = "node_id\tname\tchrom\thaplotype\tstart\tstop\tcn\n1\tx\tchr1\t1\t1\t50\t1\n"
        @test_throws ArgumentError read_profiles(IOBuffer(bad), a)   # slot incomplete
    end

    @testset "names that would corrupt a TSV are refused" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1], labels = [nothing, "a\tb"])
        @test_throws ArgumentError write_tree(IOBuffer(), t)
        res = simulate_cnas(t, toy_assembly(nchrom = 1, len = 100), CNAModel(rate = PerDivision(0.0)); seed = 1)
        @test_throws ArgumentError write_profiles(IOBuffer(), res)
        @test_throws ArgumentError write_tree(IOBuffer(),
            phylotree([nothing, 1]; labels = [nothing, "NA"]))              # would read back as missing
    end

    @testset "tree table keeps awkward values and child order" begin
        N(id, par, kids, bt, dv, mu, lb, sid) = PhyloNode(id, par, kids, bt, dv, mu, lb, sid)
        nodes = [N(1, nothing, [4, 2, 3], 0.1, nothing, nothing, "a b", nothing),
                 N(2, 1, Int[], 1 / 3, 2, 0, "say \"hi\"", 7),
                 N(3, 1, Int[], 1e-300, nothing, 5, "x'y", nothing),
                 N(4, 1, [5], nothing, 1, nothing, nothing, 9),
                 N(5, 4, Int[], 2.5e10, 3, 4, nothing, nothing)]
        t = PhyloTree(nodes)
        fields(tr) = [(n.parent, n.children, n.birthtime, n.edge_divisions,
                       n.edge_mutations, n.label, n.source_id) for n in tr.nodes]
        io = IOBuffer()
        write_tree(io, t)
        rt = read_tree(IOBuffer(take!(io)))
        @test fields(rt) == fields(t)
        @test node(rt, 1).children == [4, 2, 3]
    end

    @testset "readers validate rows" begin
        a = toy_assembly(nchrom = 1, len = 100)
        ph = "node_id\tname\tchrom\thaplotype\tstart\tstop\tcn\n"
        eh = "node_id\tname\torder\ttype\tchrom\thaplotype\tstart\tstop\tdelta\tscale\tmode\n"
        bh = "bin_index\tchrom\tstart\tstop\n"
        th = "node_id\tparent\tlabel\tsource_id\tbirthtime\tedge_divisions\tedge_mutations\n"
        r(h, body) = IOBuffer(h * body)
        @test_throws ArgumentError read_profiles(r(ph, "1\tx\tchr1\t1\t1\n"), a)        # field count
        @test_throws ArgumentError read_profiles(r(ph, "1\tx\tchrZ\t1\t1\t100\t1\n"), a)  # chromosome
        @test_throws ArgumentError read_profiles(r(ph, "1\tx\tchr1\t9\t1\t100\t1\n"), a)  # haplotype
        ev(chrom, hap, stop) = "1\tx\t1\tsegmental\t$chrom\t$hap\t1\t$stop\t1\tfocal\tNA\n"
        @test_throws ArgumentError read_events(r(eh, "1\tx\t1\n"), a)
        @test_throws ArgumentError read_events(r(eh, ev("chrZ", 1, 10)), a)
        @test_throws ArgumentError read_events(r(eh, ev("chr1", 9, 10)), a)
        @test_throws ArgumentError read_events(r(eh, ev("chr1", 1, 101)), a)
        @test length(read_events(r(eh, ev("chr1", 1, 100)), a)) == 1
        @test_throws ArgumentError read_bins(r(bh, "1\tchr1\t1\n"), a)
        @test_throws ArgumentError read_bins(r(bh, "1\tchrZ\t1\t10\n"), a)
        @test_throws ArgumentError read_bins(r(bh, "1\tchr1\t20\t10\n"), a)
        @test_throws ArgumentError read_bins(r(bh, "1\tchr1\t1\t101\n"), a)
        @test_throws ArgumentError read_tree(r(th, "1\tNA\tNA\tNA\tNA\tNA\tNA\n1\t1\tNA\tNA\tNA\tNA\tNA\n"))
        @test_throws ArgumentError read_tree(r(th, "1\tNA\tNA\tNA\tNA\tNA\tNA\n2\t7\tNA\tNA\tNA\tNA\tNA\n"))
        @test_throws ArgumentError read_tree(r(th, "1\tNA\tNA\tNA\tNA\n"))
    end

    @testset "a refused name leaves an existing file intact" begin
        a = toy_assembly(nchrom = 1, len = 100)
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1], labels = [nothing, "a\tb"])
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(2.0)); seed = 1)
        path = joinpath(mktempdir(), "x.tsv")
        write(path, "old content\n")
        @test_throws ArgumentError write_tree(path, t)
        @test_throws ArgumentError write_profiles(path, res)
        @test_throws ArgumentError write_events(path, res)
        @test read(path, String) == "old content\n"
        @test_throws ArgumentError write_events(IOBuffer(), res)
        t2 = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1], labels = [nothing, "a", "b"])
        res2 = simulate_cnas(t2, toy_assembly(nchrom = 1, len = 100), CNAModel(rate = PerDivision(1.0)); seed = 2)
        m = CNMatrix(res2, BinGrid(toy_assembly(nchrom = 1, len = 100), 50))
        @test_throws ArgumentError write_medicc2(path, m; normal_name = "a\nb")
        @test read(path, String) == "old content\n"
    end

    @testset "write_medicc2 handles a matrix with no cells" begin
        t, a, res, g = setup()
        m = CNMatrix(res, g; cells = Int[])
        buf = IOBuffer()
        write_medicc2(buf, m)
        rows = split(strip(String(take!(buf))), '\n')
        @test rows[1] == "sample_id\tchrom\tstart\tend\tcn_a\tcn_b"
        @test length(rows) - 1 == sum(length(bins_of(g, c)) for c in autosomes(a))
        @test all(startswith("diploid"), rows[2:end])
    end

    @testset "write_medicc2 refuses a cell named like the normal sample" begin
        t, a, res, g = setup()
        m = CNMatrix(res, g)
        @test_throws ArgumentError write_medicc2(IOBuffer(), m; normal_name = m.names[1])
    end
end
