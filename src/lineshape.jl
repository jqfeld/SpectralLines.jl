abstract type LineShape{T} end

# interface

(ls::LineShape)(x)= error("$(typeof(ls)) must implement the call operator.")
fwhm(ls::LineShape) = error("$(typeof(ls)) must implement a fwhm function")
hwhm(ls::LineShape) = fwhm(ls)/2

"""
    integrate(shape, a, b) -> T

Return the exact (or high-accuracy quadrature) integral ``\\int_a^b \\text{shape}(x)\\,dx``,
where `a` and `b` are positions relative to the line centre.

Use this instead of point-sampling when the grid spacing is comparable to the
linewidth (typically when `bin_width ≥ 0.3 × FWHM`).

## Methods and accuracy

| Profile | Formula | Error |
|:--------|:--------|:------|
| `Gaussian` | exact erf | < 1 ULP |
| `Lorentzian` | exact arctan | < 1 ULP |
| `PseudoVoigt` | exact (η·arctan + (1−η)·erf) | < 1 ULP |
| `Voigt` (y ≥ 0.7) | 20-pt Gauss-Hermite on Gaussian kernel | < 0.01% |
| `Voigt` (y < 0.7) | 10-pt Gauss-Legendre on [a, b] | < 0.02% for Δ ≤ 2×FWHM |

where `y = γ/(σ√2)` is the Faddeeva y-parameter stored in `Voigt.y`.

The Voigt GH branch has exact additivity (arctan telescoping); the GL branch does not —
do not use it for intervals wider than ~2 × FWHM.
"""
integrate(ls::LineShape, a, b) = error("$(typeof(ls)) must implement a integrate function")

