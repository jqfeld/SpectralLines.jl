module SpectralLines

# used in integrate(::Gaussian, a, b)
# TODO: check feasibility of pure Julia implementation to remove this dependency
import SpecialFunctions: erf 


# used in PseudoVoigt{Rocco}
# TODO: express as faddeyeva
import SpecialFunctions: erfcx

# used to calculate the weights in the quadrature for integrate(::Voigt, a, b)
using LinearAlgebra: SymTridiagonal, eigen

include("faddeyeva.jl")

include("lineshape.jl")
export fwhm, hwhm, integrate
include("lineshapes/gaussian.jl")
export Gaussian
include("lineshapes/lorentzian.jl")
export Lorentzian
include("lineshapes/voigt.jl")
export Voigt
include("lineshapes/pseudo_voigt.jl")
export PseudoVoigt, Thompson, Kielkopf, Olivero, LiuLin, Rocco



end
