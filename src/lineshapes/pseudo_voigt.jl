"""Abstract supertype for pseudo-Voigt mixing-parameter models."""
abstract type PseudoVoigtModel end

"""Thompson et al. (1987) 5th-order combined FWHM + cubic η. Default model."""
struct Thompson <: PseudoVoigtModel end

"""Kielkopf (1973) linear-quadratic combined FWHM + Thompson cubic η. Faster than Thompson."""
struct Kielkopf <: PseudoVoigtModel end

"""Olivero & Longbothum (1977) damping-factor combined FWHM + Thompson cubic η."""
struct Olivero <: PseudoVoigtModel end

"""Liu & Lin (2001) Kielkopf FWHM + cubic-in-d η (d = (fL-fG)/(fL+fG))."""
struct LiuLin <: PseudoVoigtModel end

"""Rocco & Meier (1995) erfcx-based η; η → 1 for large y. Kielkopf FWHM."""
struct Rocco <: PseudoVoigtModel end


"""
    PseudoVoigt{M<:PseudoVoigtModel, T} <: LineShape{T}

Pseudo-Voigt approximation: a linear blend ``\\eta\\,L + (1-\\eta)\\,G`` of a Lorentzian
and a Gaussian sharing a common combined HWHM.  The mixing parameter ``\\eta`` and the
combined FWHM are determined by the mixing model `M`.

Construct with `PseudoVoigt(σ, γ)` (default: `Thompson()`) or
`PseudoVoigt(σ, γ, model)` to select a specific mixing model.

Available models: [`Thompson`](@ref), [`Kielkopf`](@ref), [`Olivero`](@ref),
[`LiuLin`](@ref), [`Rocco`](@ref).

# Fields
- `sigma::T` — Gaussian width ``\\sigma``
- `gamma::T` — Lorentzian HWHM ``\\gamma``
"""
struct PseudoVoigt{M<:PseudoVoigtModel,T} <: LineShape{T}
    sigma::T
    gamma::T
end
Base.show(io::IO, pv::PseudoVoigt{M}) where {M} = print(io, "PseudoVoigt{", M, "}(σ=", pv.sigma, ", γ=", pv.gamma, ")")

PseudoVoigt(σ, γ) = PseudoVoigt(σ, γ, Thompson())
function PseudoVoigt(σ, γ, ::M) where {M<:PseudoVoigtModel}
    T = promote_type(typeof(σ), typeof(γ))
    PseudoVoigt{M,T}(T(σ), T(γ))
end
PseudoVoigt{M,T}(v::PseudoVoigt) where {M,T} = PseudoVoigt{M,T}(v.sigma, v.gamma)

function (L::PseudoVoigt{M})(x) where {M<:PseudoVoigtModel}
    α = _sqrt_2ln2 * L.sigma
    hV, η = _pseudo_voigt_params(M(), α, L.gamma)
    return η * lorentzian(x, hV) + (1 - η) * gaussian(x, hV * _inv_sqrt_2ln2)
end

@inline function fwhm(v::PseudoVoigt{M}) where {M<:PseudoVoigtModel}
    α = _sqrt_2ln2 * v.sigma
    hV, _ = _pseudo_voigt_params(M(), α, v.gamma)
    return 2 * hV
end


function integrate(pv::PseudoVoigt{M}, a, b) where {M}
    α = _sqrt_2ln2 * pv.sigma
    hV, η = _pseudo_voigt_params(M(), α, pv.gamma)
    # effective Gaussian σ for the combined width: σ_eff = hV / √(2 ln 2)
    # → 1/(σ_eff √2) = √(ln 2) / hV
    α_eff = _sqrt_ln2 / hV
    lorentz_int = (atan(b / hV) - atan(a / hV)) / π
    gauss_int = (erf(b * α_eff) - erf(a * α_eff)) / 2
    return η * lorentz_int + (1 - η) * gauss_int
end



# Implementation of the various models below


const _sqrt_2ln2 = sqrt(2 * log(2))
const _inv_sqrt_2ln2 = inv(_sqrt_2ln2)
const _sqrt_ln2 = sqrt(log(2))
const _one_minus_sqrt_pi_ln2 = 1 - sqrt(π * log(2))

@inline _thompson_eta(r) = r * evalpoly(r, (1.36603, -0.47719, 0.11116))

# _pseudo_voigt_params(model, α, γ) → (hV, η)
# α = Gaussian HWHM = √(2 ln 2)·σ;  γ = Lorentzian HWHM;  hV = combined HWHM
@inline function _pseudo_voigt_params(::Thompson, α, γ)
    hV = max(zero(α), α^5 + 2.69269 * α^4 * γ + 2.42843 * α^3 * γ^2 + 4.47163 * α^2 * γ^3 + 0.07842 * α * γ^4 + γ^5)^0.2
    return hV, _thompson_eta(γ / hV)
end

@inline function _pseudo_voigt_params(::Kielkopf, α, γ)
    hV = muladd(0.5346, γ, sqrt(muladd(0.2166, γ^2, α^2)))
    return hV, _thompson_eta(γ / hV)
end

@inline function _pseudo_voigt_params(::Olivero, α, γ)
    d = (γ - α) / (γ + α)
    hV = (1 - 0.18121 * (1 - d^2) - (0.023665 * exp(0.6 * d) + 0.00418 * exp(-1.9 * d)) * sin(π * d)) * (α + γ)
    return hV, _thompson_eta(γ / hV)
end

@inline function _pseudo_voigt_params(::LiuLin, α, γ)
    d = (γ - α) / (γ + α)
    hV = muladd(0.5346, γ, sqrt(muladd(0.2166, γ^2, α^2)))
    return hV, evalpoly(d, (0.68188, 0.61293, -0.18384, -0.11568))
end

@inline function _pseudo_voigt_params(::Rocco, α, γ)
    y = γ * _sqrt_ln2 / α      # = Faddeeva y parameter (γ / (σ√2))
    ys = min(y, oftype(y, 10))
    ex = ys * evalpoly(ys, (-0.6055, 0.0718, -0.0049, 0.000136))
    Vy = (ys + _sqrt_ln2 * exp(ex)) * erfcx(ys)
    η = y > 10 ? one(Vy) : (Vy - _sqrt_ln2) / (Vy * _one_minus_sqrt_pi_ln2)
    hV = muladd(0.5346, γ, sqrt(muladd(0.2166, γ^2, α^2)))  # Kielkopf width
    return hV, η
end


