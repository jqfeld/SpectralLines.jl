"""
    Line{T, LS<:LineShape{T}}

A single spectral line: a `shape` (any [`LineShape`](@ref)) centred at
`position` and scaled by `amplitude`. Callable as `line(x)`, evaluating
`amplitude * shape(x - position)`. Construct with `Line(position, amplitude,
shape)` — the outer constructor promotes `position`/`amplitude`/`shape`'s
numeric type parameter to a common type, so e.g. an `Int` position and a
`Float64` amplitude both end up `Float64`.

Many `Line`s are summed over a wavenumber grid with [`calculate_lines`](@ref)/
[`calculate_lines!`](@ref) to build a full spectrum.

# Fields
- `position::T` — line centre
- `amplitude::T` — scale factor
- `shape::LS` — the (possibly instrument-convolved, see [`convolve`](@ref))
  line-shape function
"""
struct Line{T,LS<:LineShape{T}}
    position::T
    amplitude::T
    shape::LS
end
Base.broadcastable(x::Line) = Ref(x)
function Base.show(io::IO, l::Line)
    print(io, "Line(pos=", l.position, ", amp=", l.amplitude, ", shape=")
    show(io, l.shape)
    print(io, ")")
end


function Line(position, amplitude, shape::LS) where {T,LS<:LineShape{T}}
    num_type = promote_type(typeof(position), typeof(amplitude), T)
    line_type = Base.typename(LS).wrapper{num_type}
    return Line{num_type,line_type}(
        position,
        amplitude,
        line_type(shape)
    )
end

@inline (l::Line{T,LS})(x) where {T,LS} = l.amplitude * l.shape(x - l.position)

integrate(l::Line, a, b) = l.amplitude * integrate(l.shape, a - l.position, b - l.position)

calculate_lines!(ret, xs, lines, method) = error("Line synthesis method $(typeof(method)) is not implemented.")
calculate_lines(xs, lines, method) = error("Line synthesis method $(typeof(method)) is not implemented.")
