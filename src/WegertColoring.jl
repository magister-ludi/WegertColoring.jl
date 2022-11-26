
module WegertColoring

using ImageCore

export zdomain, plot_phase, plot_phase!

include("phase_plot.jl")
include("colors.jl")
include("render.jl")

const description = split(
    """
WegertColoring color scheme choices
    :a => alternating black and white phase
    :b => alternating black and white modulus
    :c => phase plot with conformal polar grid
    :d => standard domain coloring
    :e => enhanced domain coloring
    :i => stripes corresponding to imaginary part
    :j => phase plot with some enhanced isochromatic lines
    :l => like i but spacing in integer fractions of 2π
    :m => enhanced phase plot with modulus jumps
    :n => like c => with brighter color for background
    :p => proper phase plot
    :q => phase plot colored in steps
    :s => conformal polar grid
    :t => polar chessboard - light gray
    :u => polar chessboard
    :v => Cartesian chessboard
    :w => cartesian chessboard - light gray
    :x => stripes corresponding to real part
    :y => stripes corresponding to imaginary part
""",
    "\n"
)

function color_schemes()
    for l in description
        println(l)
    end
end

"""
    zdomain(zll = -1.0 - im * 1.0, zur = 1.0 + im * 1.0, m = 800, n = 800)
Return a matrix of complex numbers representing the domain of a function.
`zll` and `zur` represent the lower left point and the upper right point
respectively of the domain (assuming right-hande axes).
The returned `m×n` matrix contains points `x + im * y` in the complex plane,
with `real(zll) ≤ x ≤ real(zur)` and `imag(zll) ≤ y ≤ imag(zur)`.
"""
function zdomain(zll = -1.0 - im * 1.0, zur = 1.0 + im * 1.0, m = 800, n = 800)
    a = min(real(zll), real(zur))
    b = max(real(zll), real(zur))
    # Flip y, for right-handed axes
    d = min(imag(zll), imag(zur))
    c = max(imag(zll), imag(zur))
    z = Matrix{ComplexF64}(undef, n, m)
    for col = 1:n
        zr = a + (col - 1) * (b - a) / (n - 1)
        for row = 1:m
            zi = c + (row - 1) * (d - c) / (m - 1)
            z[row, col] = zr + im * zi
        end
    end
    return z
end

end
