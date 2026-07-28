"""
    Voigt{T} <: LineShape{T}

Voigt line-shape profile: the convolution of a Gaussian (width `sigma`) and a
Lorentzian (HWHM `gamma`).  Evaluated via the real part of the Faddeeva function.

```math
V(x;\\,\\sigma,\\gamma) = \\frac{1}{\\sigma\\sqrt{2\\pi}}\\,
\\operatorname{Re}\\!\\left[w\\!\\left(\\frac{x}{\\sigma\\sqrt{2}}+i\\frac{\\gamma}{\\sigma\\sqrt{2}}\\right)\\right]
```

Normalized to unit area. Construct with `Voigt(σ, γ)`.

# Fields
- `sigma::T` — Gaussian width ``\\sigma``
- `gamma::T` — Lorentzian HWHM ``\\gamma``
"""
struct Voigt{T} <: LineShape{T}
    sigma::T
    gamma::T
    inv_sigma_sqrt2::T
    y::T
end
Voigt(s, g) = Voigt{promote_type(typeof(s), typeof(g))}(promote(s, g, inv(s * sqrt(2)), g / (s * sqrt(2)))...)
Voigt{T}(v::Voigt) where {T} = Voigt{T}(v.sigma, v.gamma, v.inv_sigma_sqrt2, v.y)


# sigma == 0 exactly makes inv_sigma_sqrt2/y = Inf, and x * Inf = NaN at x = 0
# (worse: for sigma merely *very* small but nonzero — down to ~1e-150 — the
# faddeyeva evaluation is fine, but somewhere around 1e-200 to 1e-300 it
# silently underflows to 0.0, and below ~1e-308 inv_sigma_sqrt2 itself
# overflows to Inf, giving NaN again — no value of sigma that small is ever
# physically meaningful for a real linewidth, so only the exact-zero case is
# special-cased here). gamma == 0 never crashes (it only ever appears in the
# numerator of y and inv_sigma_sqrt2 doesn't involve it at all) but is ~1e-6
# off from the exact Gaussian formula (faddeyeva's own approximation error at
# y = 0) — special-cased too since the exact formula is free.
@inline function (L::Voigt)(x)
    iszero(L.sigma) && return lorentzian(x, L.gamma)
    iszero(L.gamma) && return gaussian(x, L.sigma)
    return faddeyeva(x * L.inv_sigma_sqrt2, L.y) * L.inv_sigma_sqrt2 * inv(sqrt(pi))
end

@inline voigt(x, σ, γ) = faddeyeva(x / σ / sqrt(2), γ / σ / sqrt(2)) / σ / sqrt(2π)

# https://doi.org/10.1016/0022-4073(77)90161-3
@inline fwhm(v::Voigt{T}) where {T} = 0.5343 * v.gamma * 2 + sqrt(0.2169 * (2 * v.gamma)^2 + (2 * sqrt(2 * log(2)) * v.sigma)^2)


# 20-point Gauss-Hermite nodes/weights (physicist's convention):
# ∫_{-∞}^{∞} e^{-x²} f(x) dx ≈ Σ_k w_k f(x_k); Σ w_k = √π.
# Used by integrate(Voigt, a, b) in the Lorentzian-dominated branch (y ≥ 0.7).
const _GH20 = let
    n = 20
    J = SymTridiagonal(zeros(n), sqrt.(1:n-1) ./ sqrt(2.0))
    F = eigen(J)
    perm = sortperm(F.values)
    (nodes=F.values[perm], weights=sqrt(Float64(π)) .* F.vectors[1, perm] .^ 2)
end
const _GH20_NODES = _GH20.nodes
const _GH20_WEIGHTS = _GH20.weights
const _INV_PI_SQRTPI = inv(Float64(π) * sqrt(Float64(π)))

# 10-point Gauss-Legendre nodes/weights on (-1, 1):
# ∫_{-1}^{1} f(t) dt ≈ Σ_k w_k f(t_k); Σ w_k = 2.
# Used by integrate(Voigt, a, b) in the Gaussian-dominated branch (y < 0.7).
const _GL10 = let
    n = 10
    β = [k / sqrt(4k^2 - 1.0) for k in 1:n-1]
    J = SymTridiagonal(zeros(n), β)
    F = eigen(J)
    perm = sortperm(F.values)
    (nodes=F.values[perm], weights=2.0 .* F.vectors[1, perm] .^ 2)
end
const _GL10_NODES = _GL10.nodes
const _GL10_WEIGHTS = _GL10.weights

function integrate(v::Voigt, a, b)
    # Two-branch hybrid based on the Faddeeva y-parameter (y = γ/(σ√2)):
    #
    # y ≥ 0.7 — Lorentzian-dominated: 20-pt Gauss-Hermite on the Gaussian kernel.
    #   ∫_a^b V dx = (1/π√π) Σ_k w_k [atan((b−σ√2·x_k)/γ) − atan((a−σ√2·x_k)/γ)]
    #   Properties: exact additivity (arctan telescopes), accurate for wide intervals,
    #   error < 0.01% for all Δ when y ≥ 0.7.
    #
    # y < 0.7 — Gaussian-dominated: 10-pt Gauss-Legendre direct on [a, b].
    #   ∫_a^b V dx ≈ (b−a)/2 Σ_k w_k V(mid + half·t_k)
    #   Properties: accurate for Δ ≤ 2×FWHM, all γ/σ ratios; error < 0.012%.
    if v.y >= 0.7
        sigma_rt2 = v.sigma * sqrt(2.0)
        s = zero(promote_type(typeof(a), typeof(b), typeof(v.sigma)))
        @inbounds for k in 1:20
            t = sigma_rt2 * _GH20_NODES[k]
            s += _GH20_WEIGHTS[k] * (atan((b - t) / v.gamma) - atan((a - t) / v.gamma))
        end
        return s * _INV_PI_SQRTPI
    else
        mid = (a + b) / 2
        half = (b - a) / 2
        s = zero(promote_type(typeof(a), typeof(b), typeof(v.sigma)))
        @inbounds for k in 1:10
            s += _GL10_WEIGHTS[k] * v(mid + half * _GL10_NODES[k])
        end
        return s * half
    end
end


Base.show(io::IO, v::Voigt) = print(io, "Voigt(σ=", v.sigma, ", γ=", v.gamma, ")")
