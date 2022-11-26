
defaultjumps(pres::Integer = 20) = [exp(im * (2π * i / pres - π)) for i = 0:(pres - 1)]

"""
    plot_phase!(fz::Matrix{<:Complex}, cs::Symbol = :c, pres::Integer = 20, t::Vector{<:Complex} = defaultjumps(pres))
Return a `Matrix{RGB{N0f8} containing a plot of phases of the elements of `fz`
into the matrix `img`, using color scheme `cs`, `pres` steps around the unit
circle at locations `t`.
"""
plot_phase(
    fz::Matrix{<:Complex},
    cs::Symbol = :c,
    pres::Integer = 20,
    t::Vector{<:Complex} = defaultjumps(pres)) =
    plot_phase!(zeros(RGB{N0f8}, size(fz)), fz, cs, pres, t)

"""
    plot_phase!(img::Matrix{<:Colorant}, fz::Matrix{<:Complex}, cs::Symbol = :c, pres::Integer = 20, t::Vector{<:Complex} = defaultjumps(pres))
Plot phases of the elements of `fz` into the matrix `img`, using color
scheme `cs`, `pres` steps around the unit circle at locations `t`. Return `img`.
"""
function plot_phase!(
    img::Matrix{<:Colorant},
    fz::Matrix{<:Complex},
    cs::Symbol = :c,
    pres::Integer = 20,
    t::Vector{<:Complex} = defaultjumps(pres),
)
    @assert size(img) == size(fz)
    pres = length(t)
    return render_phase!(img, fz, cs, t, pres)
end
