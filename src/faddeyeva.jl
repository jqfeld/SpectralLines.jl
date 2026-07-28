# original implementation in https://github.comf/markmbau/Faddeyeva985.jl.git

const Θ = 1 / √π

const α₀ = 122.60793
const α₁ = 214.38239
const α₂ = 181.92853
const α₃ = 93.15558
const α₄ = 30.180142
const α₅ = 5.9126262
const α₆ = 1 / √π

const β₀ = 122.60793
const β₁ = 352.73063
const β₂ = 457.33448
const β₃ = 348.70392
const β₄ = 170.354
const β₅ = 53.992907
const β₆ = 10.479857

const γ₀ = 36183.31
const γ₁ = 3321.99
const γ₂ = 1540.787
const γ₃ = 219.031
const γ₄ = 35.7668
const γ₅ = 1.320522
const γ₆ = 1 / √π

const λ₀ = 32066.6
const λ₁ = 24322.84
const λ₂ = 9022.228
const λ₃ = 2186.181
const λ₄ = 364.2191
const λ₅ = 61.57037
const λ₆ = 1.841439

const s₀ = 38000.0
const s₁ = 256.0
const s₂ = 62.0
const s₃ = 30.0
const t₃ = 1.0e-13
const s₄ = 2.5
const t₄ = 5.0e-9
const t₅ = 0.072

#region 4: Laplace continued fractions, 4 convergents
function region4(z::Complex, x, y)::Complex
    z² = z^2
    return (Θ * (-y + im * x)) * (z² - 2.5) / (z² * (z² - 3.0) + 0.75)
end

#region 5: Humlicek's w4 (Region 4), part a
function region5a(z::Complex, x²)::Complex
    z² = z^2
    r = γ₀ + z² * (γ₁ + z² * (γ₂ + z² * (γ₃ + z² * (γ₄ + z² * (γ₅ + z² * γ₆)))))
    t = λ₀ + z² * (λ₁ + z² * (λ₂ + z² * (λ₃ + z² * (λ₄ + z² * (λ₅ + z² * (λ₆ + z²))))))
    return exp(-x²) + (im * z * r / t)
end

#region 5: Humlicek's w4 (Region 4), part b
function region5b(z::Complex)::Complex
    z² = z^2
    r = γ₀ + z² * (γ₁ + z² * (γ₂ + z² * (γ₃ + z² * (γ₄ + z² * (γ₅ + z² * γ₆)))))
    t = λ₀ + z² * (λ₁ + z² * (λ₂ + z² * (λ₃ + z² * (λ₄ + z² * (λ₅ + z² * (λ₆ + z²))))))
    return exp(-z²) + (im * z * r / t)
end

#region 6: Hui's p-6 Approximation
function region6(x, y)::Complex
    q = y - im * x
    r = α₀ + q * (α₁ + q * (α₂ + q * (α₃ + q * (α₄ + q * (α₅ + q * α₆)))))
    t = β₀ + q * (β₁ + q * (β₂ + q * (β₃ + q * (β₄ + q * (β₅ + q * (β₆ + q))))))
    return r / t
end

"""
    faddeyeva(z::Complex) -> Complex
    faddeyeva(x, y) -> Real

Evaluate the Faddeeva (scaled complex complementary error) function

```math
w(z) = e^{-z^2} \\operatorname{erfc}(-iz), \\quad z = x + iy,
```

whose real part gives the Voigt profile:
`voigt(x, σ, γ) = real(faddeyeva(x/σ/√2, γ/σ/√2)) / (σ√(2π))`.

The implementation uses a six-region piecewise rational approximation in the
``(x, y)`` plane, partitioned by ``s = x^2 + y^2``:

| Region | Condition | Method |
|:------:|-----------|--------|
| 1 | ``s \\geq 38000`` | Laplace continued fraction, 1 convergent |
| 2 | ``s \\geq 256`` | Laplace continued fraction, 2 convergents |
| 3 | ``s \\geq 62`` | Laplace continued fraction, 3 convergents |
| 4 | ``s \\geq 30``, ``y^2 \\geq 10^{-13}`` | Laplace continued fraction, 4 convergents |
| 5 | ``s > 2.5``, ``y^2 < 0.072`` | Humlicek w4 rational approximation |
| 6 | otherwise | Hui p-6 rational approximation |

The two-argument form `faddeyeva(x, y)` computes only `real(w(x + iy))` by
tracking real and imaginary parts separately through each region's polynomial
evaluation and returning only the real component. This avoids roughly one third
of the arithmetic compared to calling `real(faddeyeva(x + im*y))` and is the
form used in the hot path of [`Voigt`](@ref) evaluation.

## References

- Humlicek, J. (1982). Optimized computation of the Voigt and complex probability
  functions. *J. Quant. Spectrosc. Radiat. Transfer*, **27**(4), 437–444.
  doi:10.1016/0022-4073(82)90078-4
- Hui, A. K., Armstrong, B. H., & Wray, A. A. (1978). Rapid computation of the
  Voigt and complex error functions. *J. Quant. Spectrosc. Radiat. Transfer*,
  **19**(5), 509–516. doi:10.1016/0022-4073(78)90019-5
"""
function faddeyeva(z::Complex)::Complex

    x = real(z)
    y = imag(z)
    x² = x * x
    y² = y * y
    s = x² + y²

    #region 1: Laplace continued fractions, 1 convergent
    if s >= s₀
        return (y + im * x) * Θ / s
    end

    #region 2: Laplace continued fractions, 2 convergents
    if s >= s₁
        a = y * (0.5 + s)
        b = x * (s - 0.5)
        d = s^2 + (y² - x²) + 0.25
        return (a + im * b) * (Θ / d)
    end

    #region 3: Laplace continued fractions, 3 convergents
    if s >= s₂
        q = y² - x² + 1.5
        r = 4.0 * x² * y²
        a = y * ((q - 0.5) * q + r + x²)
        b = x * ((q - 0.5) * q + r - y²)
        d = s * (q * q + r)
        return Θ * (a + im * b) / d
    end

    #region 4: Laplace continued fractions, 4 convergents
    if s >= s₃ && y² >= t₃
        return region4(z, x, y)
    end

    #region 5: Humlicek's w4 (Region 4)
    if s > s₄ && y² < t₄
        return region5a(z, x²)
    elseif s > s₄ && y² < t₅
        return region5b(z)
    end

    #region 6: Hui's p-6 Approximation
    return region6(x, y)

end

# ── Real-only helpers for faddeyeva(x, y) ────────────────────────────────────
# The complex region functions compute both Re and Im but callers only need Re.
# These helpers track Re/Im explicitly through the polynomial evaluation and return
# only the real part, saving roughly 1/3 of the arithmetic.

function region4_real(x, y, x², y²)
    p      = x² - y²                    # Re(z²)
    q      = 2 * x * y                  # Im(z²)
    # numerator: Θ * (-y + ix) * (z² - 2.5)
    num_re = Θ * (-y * (p - 2.5) - x * q)
    num_im = Θ * ( x * (p - 2.5) - y * q)
    # denominator: z²*(z²-3) + 0.75
    den_re = p * (p - 3) - q * q + 0.75
    den_im = q * (2p - 3)
    return (num_re * den_re + num_im * den_im) / (den_re^2 + den_im^2)
end

# Horner evaluation of the γ and λ polynomials at z² = p + iq.
# Returns (Re(r), Im(r), Re(t), Im(t)) for use by region5a_real / region5b_real.
function _γλ_horner(p, q)
    r_re, r_im = γ₆, zero(p)
    r_re, r_im = r_re*p - r_im*q + γ₅, r_re*q + r_im*p
    r_re, r_im = r_re*p - r_im*q + γ₄, r_re*q + r_im*p
    r_re, r_im = r_re*p - r_im*q + γ₃, r_re*q + r_im*p
    r_re, r_im = r_re*p - r_im*q + γ₂, r_re*q + r_im*p
    r_re, r_im = r_re*p - r_im*q + γ₁, r_re*q + r_im*p
    r_re, r_im = r_re*p - r_im*q + γ₀, r_re*q + r_im*p
    t_re, t_im = λ₆ + p, q             # innermost: λ₆ + z²
    t_re, t_im = t_re*p - t_im*q + λ₅, t_re*q + t_im*p
    t_re, t_im = t_re*p - t_im*q + λ₄, t_re*q + t_im*p
    t_re, t_im = t_re*p - t_im*q + λ₃, t_re*q + t_im*p
    t_re, t_im = t_re*p - t_im*q + λ₂, t_re*q + t_im*p
    t_re, t_im = t_re*p - t_im*q + λ₁, t_re*q + t_im*p
    t_re, t_im = t_re*p - t_im*q + λ₀, t_re*q + t_im*p
    return r_re, r_im, t_re, t_im
end

# Re(exp(-x²) + im*z*r/t);  Re(im*z*r/t) = Re((-y+ix)*(r/t)) = -y*Re(r/t) - x*Im(r/t)
function region5a_real(x, y, x², y²)
    p = x² - y²
    q = 2 * x * y
    r_re, r_im, t_re, t_im = _γλ_horner(p, q)
    t² = t_re^2 + t_im^2
    return exp(-x²) + (-y * (r_re*t_re + r_im*t_im) - x * (r_im*t_re - r_re*t_im)) / t²
end

# Re(exp(-z²) + im*z*r/t);  Re(exp(-z²)) = exp(y²-x²)*cos(2xy)
function region5b_real(x, y, x², y²)
    p = x² - y²
    q = 2 * x * y
    r_re, r_im, t_re, t_im = _γλ_horner(p, q)
    t² = t_re^2 + t_im^2
    return exp(-p) * cos(q) + (-y * (r_re*t_re + r_im*t_im) - x * (r_im*t_re - r_re*t_im)) / t²
end

# Re(r/t) where r and t are Horner polynomials in q = y - ix
function region6_real(x, y)
    r_re, r_im = α₆, zero(x)
    r_re, r_im = r_re*y + r_im*x + α₅, -r_re*x + r_im*y
    r_re, r_im = r_re*y + r_im*x + α₄, -r_re*x + r_im*y
    r_re, r_im = r_re*y + r_im*x + α₃, -r_re*x + r_im*y
    r_re, r_im = r_re*y + r_im*x + α₂, -r_re*x + r_im*y
    r_re, r_im = r_re*y + r_im*x + α₁, -r_re*x + r_im*y
    r_re, r_im = r_re*y + r_im*x + α₀, -r_re*x + r_im*y
    t_re, t_im = β₆ + y, -x            # innermost: β₆ + q
    t_re, t_im = t_re*y + t_im*x + β₅, -t_re*x + t_im*y
    t_re, t_im = t_re*y + t_im*x + β₄, -t_re*x + t_im*y
    t_re, t_im = t_re*y + t_im*x + β₃, -t_re*x + t_im*y
    t_re, t_im = t_re*y + t_im*x + β₂, -t_re*x + t_im*y
    t_re, t_im = t_re*y + t_im*x + β₁, -t_re*x + t_im*y
    t_re, t_im = t_re*y + t_im*x + β₀, -t_re*x + t_im*y
    return (r_re*t_re + r_im*t_im) / (t_re^2 + t_im^2)
end

#real arguments represent z = x + im*y and return only the real part
function faddeyeva(x, y)

    x² = x * x
    y² = y * y
    s = x² + y²

    #region 1: Laplace continued fractions, 1 convergent
    if s >= s₀
        return y * Θ / s
    end

    #region 2: Laplace continued fractions, 2 convergents
    if s >= s₁
        return y * (0.5 + s) * (Θ / ((s^2 + (y² - x²)) + 0.25))
    end

    #region 3: Laplace continued fractions, 3 convergents
    if s >= s₂
        q = y² - x² + 1.5
        r = 4.0 * x² * y²
        return Θ * (y * ((q - 0.5) * q + r + x²)) / (s * (q * q + r))
    end

    if s >= s₃ && y² >= t₃
        return region4_real(x, y, x², y²)
    end

    if s > s₄ && y² < t₄
        return region5a_real(x, y, x², y²)
    elseif s > s₄ && y² < t₅
        return region5b_real(x, y, x², y²)
    end

    return region6_real(x, y)

end
