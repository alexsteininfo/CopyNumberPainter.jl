@testset "newick" begin
    @testset "branchlength must be named explicitly" begin
        @test_throws ArgumentError parse_newick("(A:1,B:2);")
        @test_throws ArgumentError parse_newick("(A:1,B:2);"; branchlength = :bogus)
    end

    @testset ":time accumulates birthtimes and nothing else" begin
        t = parse_newick("((A:1.0,B:2.0)X:0.5,C:3.0)R;"; branchlength = :time)
        @test nnodes(t) == 5
        @test node_by_label(t, "R") == treeroot(t)
        @test node(t, treeroot(t)).birthtime == 0.0
        x = node_by_label(t, "X")
        @test node(t, x).birthtime ≈ 0.5
        @test node(t, node_by_label(t, "A")).birthtime ≈ 1.5
        @test node(t, node_by_label(t, "C")).birthtime ≈ 3.0
        @test edge_time(t, node_by_label(t, "B")) ≈ 2.0
        @test all(node(t, i).edge_divisions === nothing for i in 1:nnodes(t))
        @test_throws ArgumentError n_cnas(PerDivision(1.0), t, node_by_label(t, "A"), Random.Xoshiro(1))
        @test n_cnas(PerTime(1.0), t, node_by_label(t, "A"), Random.Xoshiro(1)) isa Int
    end

    @testset ":divisions and :mutations fill their own field only" begin
        d = parse_newick("(A:3,B:5)R;"; branchlength = :divisions)
        @test node(d, node_by_label(d, "A")).edge_divisions == 3
        @test node(d, node_by_label(d, "A")).birthtime === nothing
        @test node(d, node_by_label(d, "A")).edge_mutations === nothing

        m = parse_newick("(A:3,B:5)R;"; branchlength = :mutations)
        @test node(m, node_by_label(m, "B")).edge_mutations == 5
        @test node(m, node_by_label(m, "B")).edge_divisions === nothing
        @test node(m, node_by_label(m, "B")).birthtime === nothing
    end

    @testset "structural variety" begin
        # unnamed internal nodes
        t1 = parse_newick("((A:1,B:1):1,C:2);"; branchlength = :time)
        @test length(leaves(t1)) == 3
        # multifurcation
        t2 = parse_newick("(A:1,B:1,C:1,D:1)R;"; branchlength = :time)
        @test length(childrenof(t2, treeroot(t2))) == 4
        # unary node, as a pruned lineage tree produces
        t3 = parse_newick("((A:1):2)R;"; branchlength = :time)
        @test leaves(t3) == [node_by_label(t3, "A")]
        @test depth(t3, node_by_label(t3, "A")) == 2
        # a single leaf
        t4 = parse_newick("A;"; branchlength = :divisions)
        @test nnodes(t4) == 1 && leaves(t4) == [1]
        # quoted label containing a comma and a colon
        t5 = parse_newick("('cell,1:x':1,B:1)R;"; branchlength = :time)
        @test node_by_label(t5, "cell,1:x") isa Int
        # comments are skipped
        t6 = parse_newick("(A:1[a comment],B:1)R;"; branchlength = :time)
        @test length(leaves(t6)) == 2
        # whitespace and newlines
        t7 = parse_newick("(\n  A:1 ,\n  B:1\n) R ;"; branchlength = :time)
        @test node_by_label(t7, "R") == treeroot(t7)
    end

    @testset "missing branch lengths" begin
        d = parse_newick("(A,B)R;"; branchlength = :divisions)
        @test node(d, node_by_label(d, "A")).edge_divisions === nothing
        # under :time a missing length makes that node's birthtime, and every
        # birthtime below it, unknown rather than silently zero
        t = parse_newick("((A:1,B:1)X,C:1)R;"; branchlength = :time)
        @test node(t, node_by_label(t, "X")).birthtime === nothing
        @test node(t, node_by_label(t, "A")).birthtime === nothing
        @test node(t, node_by_label(t, "C")).birthtime ≈ 1.0
    end

    @testset "malformed input is rejected" begin
        @test_throws ArgumentError parse_newick("(A:1,B:1)"; branchlength = :time)      # no semicolon
        @test_throws ArgumentError parse_newick("(A:1,B:1;"; branchlength = :time)      # unbalanced
        @test_throws ArgumentError parse_newick("(A:1,B:1));"; branchlength = :time)    # trailing paren
        @test_throws ArgumentError parse_newick("();"; branchlength = :time)            # empty branchset
        @test_throws ArgumentError parse_newick("(A:-1,B:1);"; branchlength = :time)    # negative length
        @test_throws ArgumentError parse_newick(";"; branchlength = :time)              # empty tree
        @test_throws ArgumentError parse_newick("(,A:1)R;"; branchlength = :divisions)     # leading comma
        @test_throws ArgumentError parse_newick("(A:1,)R;"; branchlength = :divisions)     # trailing comma
        @test_throws ArgumentError parse_newick("(A:1,,B:1)R;"; branchlength = :divisions) # doubled comma
    end

    @testset "fractional lengths warn under :divisions and :mutations" begin
        t = @test_logs (:warn, r"fractional") parse_newick("(A:2.4,B:1.0)R;"; branchlength = :divisions)
        @test node(t, node_by_label(t, "A")).edge_divisions == 2
    end

    @testset "round-trip in all three modes" begin
        base = phylotree([nothing, 1, 1, 2, 2];
                         birthtimes = [0.0, 1.25, 2.5, 3.0, 4.75],
                         edge_divisions = [nothing, 2, 3, 1, 4],
                         edge_mutations = [nothing, 7, 0, 5, 11],
                         labels = ["R", "X", "c", "a", "b"])
        for mode in (:time, :divisions, :mutations)
            s = newick_string(base; branchlength = mode)
            rt = parse_newick(s; branchlength = mode)
            @test nnodes(rt) == nnodes(base)
            @test [cellname(rt, i) for i in preorder(rt)] ==
                  [cellname(base, i) for i in preorder(base)]
            for i in preorder(rt)
                isroot(rt, i) && continue
                j = node_by_label(base, cellname(rt, i))
                if mode === :time
                    @test edge_time(rt, i) ≈ edge_time(base, j)
                elseif mode === :divisions
                    @test node(rt, i).edge_divisions == node(base, j).edge_divisions
                else
                    @test node(rt, i).edge_mutations == node(base, j).edge_mutations
                end
            end
        end
    end

    @testset "writing requires the field it emits" begin
        t = phylotree([nothing, 1, 1]; labels = ["R", "a", "b"])
        @test_throws ArgumentError newick_string(t; branchlength = :time)
        @test_throws ArgumentError newick_string(t; branchlength = :divisions)
    end

    @testset "labels keyword selects the emitted name" begin
        t = phylotree([nothing, 1, 1];
                      edge_divisions = [nothing, 1, 1],
                      labels = [nothing, "a", nothing],
                      source_ids = [10, 20, 30])
        @test occursin("a", newick_string(t; branchlength = :divisions, labels = :label))
        @test occursin("cell_3", newick_string(t; branchlength = :divisions, labels = :label))
        @test occursin("cell_20", newick_string(t; branchlength = :divisions, labels = :source_id))
        @test occursin("cell_2", newick_string(t; branchlength = :divisions, labels = :id))
        @test_throws ArgumentError newick_string(t; branchlength = :divisions, labels = :bogus)
    end

    @testset "file round-trip" begin
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 2, 2], labels = ["R", "a", "b"])
        path = joinpath(mktempdir(), "tree.nwk")
        write_newick(path, t; branchlength = :divisions)
        rt = read_newick(path; branchlength = :divisions)
        @test [cellname(rt, i) for i in preorder(rt)] == ["R", "a", "b"]
    end

    @testset "labels needing quotes are quoted on write" begin
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                      labels = ["R", "has,comma", "plain"])
        s = newick_string(t; branchlength = :divisions)
        @test occursin("'has,comma'", s)
        rt = parse_newick(s; branchlength = :divisions)
        @test node_by_label(rt, "has,comma") isa Int
    end

    @testset "deep trees round-trip without recursion" begin
        n = 100_000
        t = phylotree(vcat(nothing, collect(1:n - 1));
                      edge_divisions = vcat(nothing, fill(1, n - 1)))
        rt = parse_newick(newick_string(t; branchlength = :divisions); branchlength = :divisions)
        @test nnodes(rt) == n
        @test depth(rt, only(leaves(rt))) == n - 1
        @test all(node(rt, i).edge_divisions == 1 for i in 1:n if !isroot(rt, i))
    end

    @testset "node ids are assigned in preorder" begin
        t = parse_newick("((A,B)X,C)R;"; branchlength = :divisions)
        @test [cellname(t, i) for i in 1:nnodes(t)] == ["R", "X", "A", "B", "C"]
    end

    @testset "whitespace and comments around labels and lengths" begin
        t = parse_newick("( A : 1 , B [x] : 2 ) [c] R [d] ;"; branchlength = :divisions)
        @test node_by_label(t, "R") == treeroot(t)
        @test node(t, node_by_label(t, "B")).edge_divisions == 2
        @test_throws ArgumentError parse_newick("(A,B)R [open;"; branchlength = :divisions)
    end

    @testset "duplicate leaf names are refused on write; internal ones are not" begin
        dup = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                        labels = ["R", "x", "x"])
        @test_throws ArgumentError newick_string(dup; branchlength = :divisions)
        # support values as internal labels, as inferred trees carry them, still write
        sup = parse_newick("(((A:1,B:1)100:1,C:1)100:1,D:1)R;"; branchlength = :divisions)
        @test occursin("100", newick_string(sup; branchlength = :divisions))
    end


    @testset "labels with tabs or newlines are quoted and round-trip" begin
        t = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                      labels = ["R", "a\tb", "c\nd"])
        rt = parse_newick(newick_string(t; branchlength = :divisions); branchlength = :divisions)
        @test node_by_label(rt, "a\tb") isa Int
        @test node_by_label(rt, "c\nd") isa Int
    end

    @testset ".gz paths are compressed on write and decompressed on read" begin
        t = parse_newick("((A:1,B:2)C:3,D:4)R;"; branchlength = :divisions)
        path = joinpath(mktempdir(), "t.nwk.gz")
        write_newick(path, t; branchlength = :divisions)
        @test read(path)[1:2] == [0x1f, 0x8b]          # gzip magic bytes
        back = read_newick(path; branchlength = :divisions)
        @test newick_string(back; branchlength = :divisions) ==
              newick_string(t; branchlength = :divisions)
    end

    @testset "a refused write does not truncate an existing file" begin
        dup = phylotree([nothing, 1, 1]; edge_divisions = [nothing, 1, 1],
                        labels = ["R", "x", "x"])
        path = joinpath(mktempdir(), "keep.nwk")
        write(path, "keep")
        @test_throws ArgumentError write_newick(path, dup; branchlength = :divisions)
        @test read(path, String) == "keep"
    end
end
