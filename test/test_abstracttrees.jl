@testset "AbstractTrees interface" begin
    at_loaded = try
        import AbstractTrees
        true
    catch
        false
    end

    if !at_loaded
        @info "AbstractTrees.jl is not available, so the interface tests are skipped."
    else
        import AbstractTrees

        t = phylotree([nothing, 1, 1, 2, 2]; labels = ["R", "X", "c", "a", "b"])
        r = NodeRef(t)
        @test [AbstractTrees.nodevalue(n) for n in AbstractTrees.PreOrderDFS(r)] == preorder(t)
        @test [AbstractTrees.nodevalue(n) for n in AbstractTrees.Leaves(r)] == [4, 5, 3]
        @test AbstractTrees.parent(NodeRef(t, 4)).id == 2
        @test AbstractTrees.parent(r) === nothing
        @test occursin("a", sprint(AbstractTrees.print_tree, r))
        @test AbstractTrees.treeheight(r) == 2
    end
end
