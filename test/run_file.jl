# Run one test file in isolation, with the shared fixtures loaded.
#
#   julia --project=. test/run_file.jl test_newick.jl
#
# The full suite is still `Pkg.test()`; this exists for fast iteration.
using CopyNumberPainter
using Test
using Random
using Distributions

include(joinpath(@__DIR__, "fixtures.jl"))

@testset "$(ARGS[1])" begin
    include(joinpath(@__DIR__, ARGS[1]))
end
