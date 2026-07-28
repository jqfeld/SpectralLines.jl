using Test, SafeTestsets


@safetestset "Faddeyeva" begin
    include("test_faddeyeva.jl")
end

@safetestset "Lineshapes" begin
    @safetestset "Gaussian" include("test_gaussian.jl")
    @safetestset "Lorentzian" include("test_lorentzian.jl")
    @safetestset "Voigt" include("test_voigt.jl")
end

@safetestset "Code quality (Aqua.jl)" begin
    using Aqua, SpectralLines
    Aqua.test_all(SpectralLines)
end

@safetestset "Code linting (JET.jl)" begin
    include("test_jet.jl")
end
