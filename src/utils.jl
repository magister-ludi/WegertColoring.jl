
"""
    zplane(n=1200)
Return a matrix of complex numbers representing the entire complex plane
"""
function zplane(n = 1200)
    p, q, r = sphere(n)
    return stereoS2P.(p, q, r)
end

"""
    stereoP2s(z)
Stereographic projection from plane to sphere
"""
function stereoP2s(z)
    x = real(z)
    y = imag(z)

    w = x^2 + y^2 + 1

    p = 2 * x / w
    q = 2 * y / w
    r = (w - 2) / w
    return p, q, r
end

"""
    stereoS2P(x, y, z)
Stereographic projection from sphere to plane
"""
stereoS2P(x, y, z) = (x + im * y) / (1 - z)

"""
    unitcirc(n = 1024)
Return `n` points on the circumference of the unit
circle in the complex plane.
"""
unitcirc(n = 1024) = exp.(im .* [2π * i / n for i = 0:n])

"""
    sphere(n = 20)
Construct a set of 3D points on a unit sphere
with `n` faces around the circumference. Return
`(x, y, z)`, three `(n+1)×(n+1)` matrices representing
the point coordinates.
"""
function sphere(n = 10)
    x = Matrix{Float64}(undef, n + 1, n + 1)
    y = Matrix{Float64}(undef, n + 1, n + 1)
    z = Matrix{Float64}(undef, n + 1, n + 1)
    for i = 0:n
        i1 = i + 1
        θ = π * (2 * i - n) / n
        for j = 0:n
            j1 = j + 1
            ϕ = π * (2 * j - n) / 2n
            z[j1, i1] = sin(ϕ)
            x[j1, i1] = cos(ϕ) * cos(θ)
            y[j1, i1] = cos(ϕ) * sin(θ)
        end
    end
    return x, y, z
end
