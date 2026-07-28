using Test, SafeTestsets


@safetestset "Faddeyeva" begin
    include("test_faddeyeva.jl")
end

@safetestset "Lineshapes" begin
    @safetestset "Gaussian" include("test_gaussian.jl")
end

@safetestset "Code quality (Aqua.jl)" begin
    using Aqua, SpectralLines
    Aqua.test_all(SpectralLines)
end

@safetestset "Code linting (JET.jl)" begin
    include("test_jet.jl")
end
