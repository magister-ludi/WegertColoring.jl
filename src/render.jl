
function colors_a!(img, f, myangle, phaseres, _)
    ## alternating black and white stripes corresponding to phase
    black = sawfct.(myangle, 1 / phaseres, 0, 1)
    bmin, bmax = extrema(black)
    black = floor.(2 .* (black .- bmin) ./ (bmax - bmin))
    for i in eachindex(black, img, f)
        val = clamp(black[i], 0, 1)
        img[i] = construct_color(f[i], val, val, val)
    end
end

function colors_b!(img, f, _, phaseres, _)
    ## alternating black and white stripes corresponding to modulus
    black = sawfct.(log.(abs.(f)), 2π / phaseres, 0, 1)
    bmin, bmax = extrema(black)
    black = floor.(2 .* (black .- bmin) ./ (bmax - bmin))
    for i in eachindex(black, img, f)
        val = clamp(black[i], 0, 1)
        img[i] = construct_color(f[i], val, val, val)
    end
end

function colors_c!(img, f, myangle, phaseres, nphase)
    ## phase plot with conformal polar grid
    palette = hsv(600)

    blackp = sawfct.(myangle, 1 / phaseres, 0.7, 1)
    blackm = sawfct.(log.(abs.(f)), 2π / phaseres, 0.7, 1)

    #black = blackp .* blackm
    for i in eachindex(blackp, blackm, nphase, img, f)
        black = blackp[i] * blackm[i]
        if nphase[i] > length(palette)
            img[i] = RGB(1, 1, 1)
        else
            c = palette[nphase[i]]
            img[i] = construct_color(f[i], black * red(c), black * green(c), black * blue(c))
        end
    end
end
