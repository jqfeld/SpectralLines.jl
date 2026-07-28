"""
    Lorentzian{T} <: LineShape{T}

Cauchy (Lorentzian) line-shape profile parameterized by its half-width at half-maximum.

```math
L(x;\\,\\gamma) = \\frac{1}{\\pi}\\frac{\\gamma}{x^2 + \\gamma^2}
```

Normalized to unit area. Construct with `Lorentzian(γ)`.

# Fields
- `hwhm::T` — half-width at half-maximum ``\\gamma``
"""
struct Lorentzian{T} <: LineShape{T}
    hwhm::T
end
Lorentzian{T}(l::Lorentzian) where {T} = Lorentzian{T}(l.hwhm) # copy constructor

@inline lorentzian(x, γ) = γ / ((abs2(γ) + abs2(x)) * π)
@inline (L::Lorentzian{T})(x) where {T} = lorentzian(x, L.hwhm)
@inline fwhm(l::Lorentzian{T}) where {T} = 2 * l.hwhm
integrate(l::Lorentzian, a, b) = (atan(b / l.hwhm) - atan(a / l.hwhm)) / π

Base.show(io::IO, l::Lorentzian) = print(io, "Lorentzian(γ=", l.hwhm, ")")




