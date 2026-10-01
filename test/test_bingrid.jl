@testset "bingrid" begin
    S = CopyNumberEvolution.Segment

    @testset "grid layout" begin
        a = toy_assembly(nchrom = 2, len = 1000)
        g = BinGrid(a, 250)
        @test nbins(g) == 8
        @test bins_of(g, 1) == 1:4
        @test bins_of(g, 2) == 5:8
        @test g.bins[1] == CopyNumberEvolution.Bin(1, 1, 250)
        @test g.bins[4] == CopyNumberEvolution.Bin(1, 751, 1000)
        @test_throws ArgumentError BinGrid(a, 0)
    end

    @testset "the last bin of a chromosome is short when the length does not divide" begin
        a = toy_assembly(nchrom = 1, len = 1000)
        g = BinGrid(a, 300)
        @test nbins(g) == 4
        @test g.bins[4] == CopyNumberEvolution.Bin(1, 901, 1000)
        @test length(bins_of(g, 1)) == 4
    end

    @testset "hg38 bin counts" begin
        g = BinGrid(hg38(:female), 500_000)
        @test length(bins_of(g, 1)) == cld(248956422, 500_000)   # 498
        y = chromindex(hg38(:female), "chrY")
        @test isempty(bins_of(g, y))                              # no slots, no bins
        gm = BinGrid(hg38(:male), 500_000)
        @test length(bins_of(gm, chromindex(hg38(:male), "chrY"))) == cld(57227415, 500_000)
    end

    @testset "a bin-aligned segmentation projects exactly" begin
        a = toy_assembly(nchrom = 1, len = 1000)
        g = BinGrid(a, 250)
        segs = [S(1, 250, 1), S(251, 500, 3), S(501, 1000, 0)]
        @test project(segs, g, 1, LengthWeightedMajority()) == [1, 3, 0, 0]
        @test project(segs, g, 1, AreaWeightedMean()) == [1, 3, 0, 0]
    end

    @testset "a chromosome-length segment fills every bin" begin
        a = toy_assembly(nchrom = 1, len = 1000)
        g = BinGrid(a, 250)
        @test project([S(1, 1000, 4)], g, 1, LengthWeightedMajority()) == [4, 4, 4, 4]
    end

    @testset "a straddling bin follows the documented rule" begin
        a = toy_assembly(nchrom = 1, len = 400)
        g = BinGrid(a, 100)
        # bin 2 is 101:200; the breakpoint at 161 gives 60 bp at cn 1 and 40 bp at cn 4
        segs = [S(1, 160, 1), S(161, 400, 4)]
        @test project(segs, g, 1, LengthWeightedMajority()) == [1, 1, 4, 4]
        # area-weighted mean over bin 2 is (60*1 + 40*4)/100 = 2.2 -> 2
        @test project(segs, g, 1, AreaWeightedMean()) == [1, 2, 4, 4]
    end

    @testset "majority ties resolve to the lower copy number" begin
        a = toy_assembly(nchrom = 1, len = 200)
        g = BinGrid(a, 100)
        segs = [S(1, 50, 5), S(51, 200, 2)]     # bin 1 is 50 bp of cn 5 and 50 bp of cn 2
        @test project(segs, g, 1, LengthWeightedMajority())[1] == 2
    end

    @testset "area-weighted mean rounds halves upward" begin
        a = toy_assembly(nchrom = 1, len = 100)
        g = BinGrid(a, 100)
        segs = [S(1, 50, 1), S(51, 100, 2)]     # mean exactly 1.5
        @test project(segs, g, 1, AreaWeightedMean()) == [2]
    end

    @testset "profile projection returns total and per-haplotype tracks" begin
        a = toy_assembly(nchrom = 1, len = 400)
        g = BinGrid(a, 100)
        p = diploid(a)
        apply!(p, SegmentalCNA(1, 1, 1, 200, 1, :focal))
        tot, alleles = project(p, g)
        @test length(alleles) == 2
        @test alleles[1] == [2, 2, 1, 1]
        @test alleles[2] == [1, 1, 1, 1]
        @test tot == [3, 3, 2, 2]
        @test tot == alleles[1] .+ alleles[2]
    end

    @testset "hemizygous chromosomes give zeros on the absent haplotype" begin
        a = toy_sex_assembly(:male, len = 400)
        g = BinGrid(a, 100)
        p = diploid(a)
        tot, alleles = project(p, g)
        x = bins_of(g, chromindex(a, "chrX"))
        @test all(alleles[1][i] == 1 for i in x)
        @test all(alleles[2][i] == 0 for i in x)
        @test all(tot[i] == 1 for i in x)
    end

    @testset "CNMatrix over an evolution result" begin
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                      labels = [nothing, "left", nothing])
        a = toy_assembly(nchrom = 2, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(1.0)); seed = 31)
        g = BinGrid(a, 100)
        m = CNMatrix(res, g)
        @test ncells(m) == 2
        @test m.cells == leaves(t)
        @test m.names == ["left", "cell_3"]
        @test size(m.total) == (2, nbins(g))
        @test length(m.allele) == 2
        @test m.total == m.allele[1] .+ m.allele[2]
        @test all(>=(0), m.total)
        @test max_cn(m) >= 0
    end

    @testset "CNMatrix can include internal nodes as a truth set" begin
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(1.0)); seed = 32)
        g = BinGrid(a, 100)
        m = CNMatrix(res, g; cells = collect(1:nnodes(t)))
        @test ncells(m) == 3
        @test m.names == ["cell_1", "cell_2", "cell_3"]
    end

    @testset "allele = false omits the per-haplotype tracks" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(1.0)); seed = 33)
        m = CNMatrix(res, BinGrid(a, 100); allele = false)
        @test m.allele === nothing
        @test size(m.total, 1) == 1
    end

    @testset "grid and result must agree on the assembly" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 1, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0)); seed = 34)
        other = BinGrid(toy_assembly(nchrom = 2, len = 400), 100)
        @test_throws ArgumentError CNMatrix(res, other)
    end

    @testset "matrix total is the sum of haplotype projections, not a reprojected total" begin
        # Projection is non-linear, so a length-weighted majority of a sum is not the
        # sum of the majorities. This fixture makes the two definitions disagree:
        #   hap1: 1:60 at cn 3, 61:100 at cn 0   -> majority 3 (60bp beats 40bp)
        #   hap2: 1:40 at cn 0, 41:100 at cn 2   -> majority 2 (60bp beats 40bp)
        #   additive total                        = 3 + 2 = 5
        #   summed segmentation: (1:40,3) (41:60,5) (61:100,2), lengths 40/20/40,
        #     so the majority ties at 40bp between cn 3 and cn 2 and resolves low -> 2
        a = toy_assembly(nchrom = 1, len = 100)
        g = BinGrid(a, 100)                      # exactly one bin, 1:100
        p = diploid(a)
        apply!(p, SegmentalCNA(1, 1, 1, 60, 2, :focal))       # hap1 1:60 -> cn 3
        apply!(p, SegmentalCNA(1, 1, 61, 100, -1, :focal))    # hap1 61:100 -> cn 0
        apply!(p, SegmentalCNA(1, 2, 1, 40, -1, :focal))      # hap2 1:40 -> cn 0
        apply!(p, SegmentalCNA(1, 2, 41, 100, 1, :focal))     # hap2 41:100 -> cn 2

        tot, alleles = project(p, g)
        @test alleles[1] == [3]
        @test alleles[2] == [2]
        @test tot == [5]
        @test tot == alleles[1] .+ alleles[2]

        # and the rejected definition really does give a different answer here,
        # so this test would fail if project ever reprojected the summed segmentation
        @test project(total_cn(p, 1), g, 1, LengthWeightedMajority()) == [2]
    end

    @testset "CNMatrix refuses duplicate row names" begin
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                      labels = [nothing, "cell_3", nothing])     # node 3 is also cell_3
        a = toy_assembly(nchrom = 1, len = 400)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0)); seed = 37)
        err = try CNMatrix(res, BinGrid(a, 100)); nothing catch e; e end
        @test err isa ArgumentError && occursin("cell_3", err.msg)
        @test_throws ArgumentError CNMatrix(res, BinGrid(a, 100); cells = [2, 2])
    end

    @testset "a mask drops bins above the masked fraction" begin
        a = toy_assembly(nchrom = 1, len = 1000)
        g = BinGrid(a, 100; mask = ["chr1" => 151:260])
        # 101:200 is 50 % masked (kept: not above 0.5); 201:300 is 60 % (dropped)
        @test [(b.start, b.stop) for b in g.bins] == [(s, s + 99) for s in 1:100:901 if s != 201]
        @test g.nmasked == 1
        @test bins_of(g, 1) == 1:9
        strict = BinGrid(a, 100; mask = ["chr1" => 151:260], max_masked_fraction = 0.0)
        @test !any(b -> b.start in (101, 201), strict.bins)
        @test_throws ArgumentError BinGrid(a, 100; max_masked_fraction = 1.5)
        @test BinGrid(a, 100).nmasked == 0                       # no mask, nothing dropped
    end

    @testset "overlapping and out-of-range mask intervals are merged and clipped" begin
        a = toy_assembly(nchrom = 1, len = 1000)
        g = BinGrid(a, 100; mask = ["chr1" => 250:320, "chr1" => 300:420, "chr1" => 990:5000])
        @test g.mask == ["chr1" => 250:420, "chr1" => 990:1000]
    end

    @testset "centromere_mask removes centromeric bins" begin
        a = toy_assembly(nchrom = 2, len = 1000)                # centromeres 401:600
        @test centromere_mask(a) == ["chr1" => 401:600, "chr2" => 401:600]
        g = BinGrid(a, 100; mask = centromere_mask(a))
        @test nbins(g) == 16
        @test !any(b -> 401 <= b.start <= 600, g.bins)
    end

    @testset "read_bed_mask converts BED to 1-based inclusive" begin
        bed = "# comment\ntrack name=x\nbrowser position chr1\n" *
              "chr1\t150\t260\tblack\r\nchrM 0 100\n\n"         # CRLF and space-separated
        m = read_bed_mask(IOBuffer(bed))
        @test m == ["chr1" => 151:260, "chrM" => 1:100]
        a = toy_assembly(nchrom = 1, len = 1000)
        g = @test_logs (:warn, r"not in assembly") BinGrid(a, 100; mask = m)
        @test g.nmasked == 1
        @test_throws ArgumentError read_bed_mask(IOBuffer("chr1\t10\n"))
        @test_throws ArgumentError read_bed_mask(IOBuffer("chr1\t20\t10\n"))
        @test_throws ArgumentError read_bed_mask(IOBuffer("chr1\tx\t10\n"))
        path = joinpath(mktempdir(), "m.bed")
        write(path, "chr1\t0\t10\n")
        @test read_bed_mask(path) == ["chr1" => 1:10]
        gz = joinpath(mktempdir(), "m.bed.gz")
        CopyNumberEvolution._with_io(io -> print(io, "chr1\t0\t10\nchr2\t5\t9\n"), gz)
        @test read_bed_mask(gz) == ["chr1" => 1:10, "chr2" => 6:9]
        # only a whole first field of `track` / `browser` marks a header
        @test read_bed_mask(IOBuffer("trackA\t0\t10\nbrowserB\t0\t5\ntrack\tx\n")) ==
              ["trackA" => 1:10, "browserB" => 1:5]
        @test_throws ArgumentError read_bed_mask(IOBuffer("chr1\t10\t10\n"))   # zero-length
    end

    @testset "masked grids project and export" begin
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        a = toy_assembly(nchrom = 2, len = 1000)
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(2.0)); seed = 35)
        g = BinGrid(a, 100; mask = centromere_mask(a))
        m = CNMatrix(res, g)
        @test size(m.total) == (1, 16)
        buf = IOBuffer()
        write_medicc2(buf, m)
        @test count(==('\n'), String(take!(buf))) == 1 + 2 * 16   # header, normal, one cell
        @test occursin("masked", sprint(show, g))
    end

    @testset "major_minor sorts alleles per bin and preserves totals" begin
        a = toy_assembly(nchrom = 1, len = 400)
        start = diploid(a)
        apply!(start, SegmentalCNA(1, 1, 1, 200, 2, :focal))       # hap1 1:200 -> 3
        apply!(start, SegmentalCNA(1, 2, 201, 400, 1, :focal))     # hap2 201:400 -> 2
        t = phylotree([nothing, 1]; edge_divisions = [nothing, 1])
        res = simulate_cnas(t, a, CNAModel(rate = PerDivision(0.0), initial = Given(start)); seed = 36)
        m = CNMatrix(res, BinGrid(a, 100))
        @test m.phasing === :haplotype
        @test m.allele[1] == [3 3 1 1] && m.allele[2] == [1 1 2 2]
        u = major_minor(m)
        @test u.phasing === :major_minor
        @test u.allele[1] == [3 3 2 2] && u.allele[2] == [1 1 1 1]
        @test u.total == m.total
        @test major_minor(u).allele == u.allele                     # idempotent
        @test m.allele[1] == [3 3 1 1]                              # input untouched
        @test occursin("major/minor", sprint(show, u))
        @test_throws ArgumentError major_minor(CNMatrix(res, BinGrid(a, 100); allele = false))
    end

    @testset "the sweep projection equals per-bin search" begin
        rng = Random.Xoshiro(38)
        a = toy_assembly(nchrom = 1, len = 5000)
        for _ in 1:50
            p = diploid(a)
            for _ in 1:15
                s = rand(rng, 1:5000); e = min(5000, s + rand(rng, 0:800))
                apply!(p, SegmentalCNA(1, 1, s, e, rand(rng, (-1, 1, 2)), :focal))
            end
            segs = slot_segments(p, 1, 1)
            for size in (97, 250, 1000), rule in (LengthWeightedMajority(), AreaWeightedMean())
                g = BinGrid(a, size; mask = ["chr1" => 1200:1900])
                slow = [CopyNumberEvolution._bin_value(segs, CopyNumberEvolution.segment_index(segs, b.start), b, rule)
                        for b in g.bins]
                @test project(segs, g, 1, rule) == slow
            end
        end
    end

    @testset "CNMatrix works on a lean run" begin
        t = phylotree([nothing, 1, 1, 2, 2]; edge_divisions = [nothing, 1, 1, 1, 1])
        a = toy_assembly(nchrom = 2, len = 1000)
        m = CNAModel(rate = PerDivision(2.0))
        full = simulate_cnas(t, a, m; seed = 39)
        lean = simulate_cnas(t, a, m; seed = 39, retain_internal = false)
        g = BinGrid(a, 100)
        ids = collect(1:nnodes(t))
        @test CNMatrix(lean, g; cells = ids).allele == CNMatrix(full, g; cells = ids).allele
    end
end
