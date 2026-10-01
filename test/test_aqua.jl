using Aqua

@testset "Aqua" begin
    # Ambiguities are checked for this package only; Distributions and StatsBase
    # carry their own, which are not ours to fix.
    Aqua.test_all(CopyNumberPainter; ambiguities = (recursive = false,))
end
