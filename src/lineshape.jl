abstract type LineShape{T} end

# interface

(ls::LineShape)(x)= error("$(typeof(ls)) must implement the call operator.")
fwhm(ls::LineShape) = error("$(typeof(ls)) must implement a fwhm function")
hwhm(ls::LineShape) = fwhm(ls)/2

integrate(ls::LineShape, a, b) = error("$(typeof(ls)) must implement a integrate function")

