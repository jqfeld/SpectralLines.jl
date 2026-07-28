module SpectralLines

# used in integrate(::Gaussian, a, b)
# TODO: check feasibility of pure Julia implementation to remove this dependency
import SpecialFunctions: erf 

include("faddeyeva.jl")

include("lineshape.jl")
export fwhm, hwhm, integrate
include("lineshapes/gaussian.jl")
export Gaussian


end
