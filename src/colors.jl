
#=
Methods for mapping a complex valued matrix to a matix
of color types.
=#

"""
    hsv([::Type{T} where {T <: Colorant} = RGB{Float64},] n::Integer)
Return a `Vector{T}` containing `n` values with equally space hue values
spanning `0 ≤ hue ≤ 360°` with maximum saturation and value.
"""
hsv(n) = hsv(RGB{Float64}, n)

hsv(::Type{T}, n) where {T <: Colorant} = T.([HSV(i * 360 / n, 1, 1) for i = 0:(n - 1)])

function render_phase!(img, f, cs = :p, t = nothing, phaseres = 20)
    if t === nothing || isempty(t)
        t = exp.([2π * im * i / 20 - π for i = 0:19])
    end
    m, n = size(f)

    ## choose basic color map

    # stepwise coloring or continuous coloring?
    colmap = cs == :q ? hsv(phaseres) : hsv(600)

    # number of colors in palette encoding phase
    p = length(colmap)
    # encoding phase
    myangle = (angle.(-f) .+ π) ./ 2π
    nphase = stepfct.(myangle, p)

    if cs == :a
        colors_a!(img, f, myangle, phaseres, nphase)
        #=
        ## alternating black and white stripes corresponding to phase
        black = sawfct.(myangle, 1 / phaseres, 0, 1)
        bmin, bmax = extrema(black)
        black = floor.(2 .* (black .- bmin) ./ (bmax - bmin))
        for i = 1:length(black)
            val = clamp(black[i], 0, 1)
            img[i] = construct_color(f[i], val, val, val)
        end
        =#
    elseif cs == :b
        colors_b!(img, f, myangle, phaseres, nphase)
        #=
        ## alternating black and white stripes corresponding to modulus

        black = sawfct.(log.(abs.(f)), 2π / phaseres, 0, 1)

        bmin, bmax = extrema(black)
        black = floor.(2 .* (black .- bmin) ./ (bmax - bmin))
        for i = 1:length(black)
            val = clamp(black[i], 0, 1)
            img[i] = construct_color(f[i], val, val, val)
        end
        =#
    elseif cs == :c
        colors_c!(img, f, myangle, phaseres, nphase)
        #=
        ## phase plot with conformal polar grid
        palette = hsv(600)

        blackp = sawfct.(myangle, 1 / phaseres, 0.7, 1)
        blackm = sawfct.(log.(abs.(f)), 2π / phaseres, 0.7, 1)

        #black = blackp .* blackm
        for i = 1:length(blackp)
            black = blackp[i] * blackm[i]
            if nphase[i] > length(palette)
                img[i] = RGB(1, 1, 1)
            else
                c = palette[nphase[i]]
                img[i] = construct_color(f[i], black * red(c), black * green(c), black * blue(c))
            end
        end
        =#
    elseif cs == :d || cs == :e
        ##  proper and enhanced domain coloring
        palette = hsv(600)
        if cs == :d
            for i = 1:length(img)
                img[i] = palette[nphase[i]]
            end
        elseif cs == :e
            black = sawfct.(log.(abs.(f)), 2π / phaseres, 0.7, 1)
            for i = 1:length(black)
                c = palette[nphase[i]]
                img[i] = construct_color(
                    f[i],
                    black[i] * red(c),
                    black[i] * green(c),
                    black[i] * blue(c),
                )
            end
        end

        # modified domain coloring
        bright = (3.2 / π) * atan.(abs.(f)) .- 0.8
        img .= brighten.(img, bright)
    elseif cs == :g
        base = π
        scale = 1 / 2π
        for i = 1:length(img)
            val =  scale * (base + angle(f[i]))
            img[i] = construct_color(1, val, val, val)
        end
    elseif cs == :i || cs == :l
        ## alternating black and white stripes corresponding to imaginary part

        impart = imag.(f)

        fmin, fmax = extrema(impart)
        imres = (fmax - fmin) / phaseres

        if cs == :l
            black = sawfct.(impart, 2 * π / phaseres, 0, 1)
        else
            black = sawfct.(impart, imres, 0, 1)
        end

        bmin, bmax = extrema(black)
        black = floor.(2 * (black .- bmin) ./ (bmax - bmin))
        for i = 1:length(black)
            val = clamp(black[i], 0, 1)
            img[i] = construct_color(f[i], val, val, val)
        end
    elseif cs == :j
        ## phase plot with specified jumps at points t

        palette = hsv(600)

        t = exp.(im .* sort(angle.(t)))
        black = sawfctt(f, t)

        for i = 1:length(black)
            c = palette[nphase[i]]
            img[i] =
                construct_color(f[i], black[i] * red(c), black[i] * green(c), black[i] * blue(c))
        end
    elseif cs == :m
        ## contour lines and conformal grid

        palette = hsv(600)
        black = sawfct.(log.(abs.(f)), 2π / phaseres, 0.7, 1)
        for i = 1:length(black)
            c = palette[nphase[i]]
            img[i] =
                construct_color(f[i], black[i] * red(c), black[i] * green(c), black[i] * blue(c))
        end
    elseif cs == :n
        ## phase plot with conformal polar grid for background
        palette = hsv(600)

        blackp = sawfct.(myangle, 1 / phaseres, 0.7, 1)
        blackm = sawfct.(log.(abs.(f)), 2π / phaseres, 0.7, 1)

        black = blackp .* blackm
        for i = 1:length(black)
            c = palette[nphase[i]]
            img[i] =
                construct_color(f[i], black[i] * red(c), black[i] * green(c), black[i] * blue(c))
        end
        img .= brighten.(img, 0.6)
    elseif cs == :p
        ## proper phase portrait
        palette = hsv(600)
        for i = 1:length(img)
            img[i] = palette[nphase[i]]
        end
    elseif cs == :q
        ##  phase plot colored in steps

        palette = hsv(phaseres)
        for i = 1:length(img)
            img[i] = palette[nphase[i]]
        end
    elseif cs == :s
        ## special conformal grid

        phase = sawfct.(myangle, 1 / phaseres, 0, 1)
        modul = sawfct.(log.(abs.(f)), 2π / phaseres, 0, 1)
        for i = 1:length(img)
            img[i] = construct_color(
                f[i],
                phase[i] * (1 - modul[i]),
                modul[i] * (1 - phase[i]),
                phase[i] * modul[i],
            )
        end
    elseif cs == :t
        ## polar chessboard - light gray

        blackp = rectfct.(myangle, 1 / phaseres)
        blackm = rectfct.(log.(abs.(f)), 2π / phaseres)

        black = mod.(blackp .+ blackm, 2)
        for i = 1:length(img)
            val = 0.85 + 0.15 * black[i]
            img[i] = construct_color(f[i], val, val, val)
        end
    elseif cs == :u
        ## polar chessboard
        blackp = rectfct.(myangle, 1 / phaseres)
        blackm = rectfct.(log.(abs.(f)), 2π / phaseres)

        black = mod.(blackp + blackm, 2)
        for i = 1:length(img)
            img[i] = construct_color(f[i], black[i], black[i], black[i])
        end
    elseif cs == :v
        ## cartesian chessboard
        sc = 10
        blackx = rectfct.(real.(f), sc / phaseres)
        blacky = rectfct.(imag.(f), sc / phaseres)
        white = mod.(blackx .+ blacky, 2)
        for i = 1:length(img)
            img[i] = construct_color(f[i], white[i], white[i], white[i])
        end
    elseif cs == :w
        ## cartesian chessboard, light gray
        sc = 10
        blackx = rectfct.(real.(f), sc / phaseres)
        blacky = rectfct.(imag.(f), sc / phaseres)
        black = mod.(blackx .+ blacky, 2)
        for i = 1:length(img)
            val = 0.85 + 0.15 * black[i]
            img[i] = construct_color(f[i], val, val, val)
        end
    elseif cs == :x
        ## alternating black and white stripes corresponding to real part
        reres = 10 / phaseres
        black = sawfct.(real.(f), reres, 0, 1)
        bmin, bmax = extrema(black)
        black = floor.(2 .* (black .- bmin) ./ (bmax - bmin))
        for i = 1:length(img)
            val = clamp(black[i], 0, 1)
            img[i] = construct_color(f[i], val, val, val)
        end
    elseif cs == :y
        ## alternating black and white stripes corresponding to imaginary part
        imres = 10 / phaseres
        if cs == :l
            black = sawfct.(imag.(f), 2π / pres, 0, 1)
        else
            black = sawfct.(imag.(f), imres, 0, 1)
        end
        bmin, bmax = extrema(black)
        black = floor.(2 .* (black .- bmin) ./ (bmax - bmin))
        for i = 1:length(img)
            val = clamp(black[i], 0, 1)
            img[i] = construct_color(f[i], val, val, val)
        end
    else
        println("Unknown colour scheme '$cs'")
        exit(1)
    end
    return img
end

function construct_color(f, r, g, b)
    ## coloring zeros black, poles white ######################################
    if iszero(f)
        RGB(0, 0, 0)
    elseif isnan(f) || isinf(f)
        RGB(1, 1, 1)
    else
        RGB(r, g, b)
    end
end

# saw tooth function on R with period dx onto [a,b]
function sawfct(x, dx, a, b)
    x = x / dx - floor(x / dx)
    return a + (b - a) * x
end

# integer step function of x with period 1 such that [0,1] --> [1,nmax]
function stepfct(x, nmax)
    if isinf(x) || isnan(x)
        return typemax(Int)
    else
        x = x - floor(x)
        y = floor(Int, nmax * x) + 1
        return y + (y == 0)
    end
end

function sawfctx(x, xk)
    # saw tooth function on [0,1] with values in [0,1] and jumps at points xk

    xk = sort(xk)
    kk = length(xk)
    xk = vcat(xk[kk] - 1, xk, xk[1] + 1)
    y = zeros(size(x))

    for j = 1:(kk + 1)
        for l = 1:length(x)
            if xk[j] ≤ x[l] < xk[j + 1]
                y[l] += (x[l] - xk[j]) / (xk[j + 1] - xk[j])
            end
        end
    end

    return y .- floor.(y)
end

## auxiliary function: rectfct rectangular impulses of length dx
rectfct(x, dx) = mod(floor(x / dx), 2)

## auxiliary function: sawfctt
function sawfctt(f, t, min = 0.8, max = 1.0)
    # sawtooth function on unit circle with jumps min to max at points t
    #
    # f - function to be represented
    # t - vector containing jump points on unit circle
    # min - min value of brightness
    # max - max value of brightness

    # eliminate multiple values of t
    tau = unique(sort(angle.(t)))
    t = exp.(im .* tau)

    ang = (-angle.(f) .+ π) ./ 2π
    angjmp = (angle.(conj.(t)) .+ π) ./ 2π
    saw = sawfctx(ang, angjmp)
    return min .+ (max - min) .* saw
end

function brighten(rgb, bright)
    # modification of color scheme
    # bright - scalar value or field
    # between 0 and 1 for brightening
    # between -1 and 0 for darkening
    @assert -1 ≤ bright ≤ 1
    r, g, b = red(rgb), green(rgb), blue(rgb)
    if bright > 0
        return RGB((1 - bright) * r + bright, (1 - bright) * g + bright, (1 - bright) * b + bright)
    else
        return RGB((1 + bright) * r, (1 + bright) * g, (1 + bright) * b)
    end
end
