"""
    Gaussian{T} <: LineShape{T}

Gaussian line-shape profile parameterized by its standard deviation `sigma`.

```math
G(x;\\,\\sigma) = \\frac{1}{\\sigma\\sqrt{2\\pi}}\\exp\\left(-\\frac{x^2}{2\\sigma^2}\\right)
```

Normalized to unit area. Construct with `Gaussian(σ)`.

# Fields
- `sigma::T` — standard deviation ``\\sigma``
"""
struct Gaussian{T} <: LineShape{T}
    sigma::T
    norm::T
    alpha::T
end
Gaussian(σ) = Gaussian{promote_type(typeof(σ), typeof(inv(sqrt(2π) * σ)))}(σ, inv(sqrt(2π) * σ), inv(σ * sqrt(2)))
Gaussian{T}(σ) where {T} = Gaussian{T}(σ, inv(sqrt(2π) * σ), inv(σ * sqrt(2)))
Gaussian{T}(g::Gaussian) where {T} = Gaussian{T}(g.sigma, g.norm, g.alpha)

# interface 
@inline (L::Gaussian)(x) = L.norm * exp(-abs2(x * L.alpha))
@inline fwhm(g::Gaussian{T}) where {T} = 2 * sqrt(2 * log(2)) * g.sigma
integrate(g::Gaussian, a, b) = (erf(b * g.alpha) - erf(a * g.alpha)) / 2

# plotting interface
Base.show(io::IO, g::Gaussian) = print(io, "Gaussian(σ=", g.sigma, ")")


# Not used for the line profile calculation itself, only a convenience function
@inline gaussian(x, σ) = 1 / (sqrt(2π) * σ) * exp(-abs2(x / (σ * sqrt(2))))
