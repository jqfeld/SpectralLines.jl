using Test, SafeTestsets


@safetestset "Faddeyeva" begin
    include("test_faddeyeva.jl")
end

@safetestset "Code quaity (Aqua.jl)" begin
    using Aqua, SpectralLines
    Aqua.test_all(SpectralLines)
end

@safetestset "Code linting (JET.jl)" begin
    include("test_jet.jl")
end
