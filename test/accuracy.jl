using Random
using StatsBase: std

@testset "accuracy.jl" begin
    # `μ` is the local mean, sign included: a sky-subtracted image has negative ones,
    # which `sqrt(μ1μ2)` returned as positive.
    Np, wid = 3, 5
    img = fill(-2.0, 20, 20)
    bimage = zeros(16, 16)
    bism = zeros(16, 16, 2Np - 1, Np)
    cov_avg!(bimage, similar(img), bism, img; Np, widx = wid, widy = wid, ftype = 64)
    μ = zeros(Np^2)
    build_cov!(zeros(Np^2, Np^2), μ, 8, 8, bimage, bism, Np, wid, wid)
    @test all(≈(-2.0), μ)
    bism_sym = zeros(16, 16, 2Np - 1, 2Np - 1)
    cov_avg_sym!(bimage, similar(img), bism_sym, img; Np, widx = wid, widy = wid, ftype = 64)
    build_cov_sym!(zeros(Np^2, Np^2), μ, 8, 8, bimage, bism_sym, Np, wid, wid)
    @test all(≈(-2.0), μ)

    # Float32 infill on a bright background must match Float64.  The covariance is
    # E[xy] - E[x]E[y], which cancels catastrophically in Float32 when the mean is
    # large against the fluctuations (here mean²/variance ~ 2e7).
    rng = MersenneTwister(1)
    n = 80
    field = [3 * sin(i / 5) * cos(j / 7) for i in 1:n, j in 1:n] .+ 0.3 .* randn(rng, n, n)
    img = 1.0e4 .+ field
    mask = falses(n, n)
    mask[39:41, 39:41] .= true
    out32 = proc_discrete([40], [40], Float32.(img), copy(mask); Np = 9, widx = 25, ftype = 32)
    out64 = proc_discrete([40], [40], copy(img), copy(mask); Np = 9, widx = 25, ftype = 64)
    @test maximum(abs.(out32[mask] .- out64[mask])) < 0.05

    # The caller's image is left untouched
    img32 = Float32.(img)
    for f in (im -> proc_discrete([40], [40], im, copy(mask); Np = 9, widx = 25),
              im -> proc_continuous(im, copy(mask); Np = 9, widx = 25))
        im = copy(img32)
        f(im)
        @test im == img32
    end

    # On pure white noise (σ = 1) with 30% of pixels masked, draws for the masked pixels
    # must scatter about their mean by σ.  That holds only if the noise injected into
    # masked pixels for covariance training matches the image's own white-noise level,
    # i.e. is estimated from unmasked adjacent pairs and scaled by 1/√2.
    rng = MersenneTwister(1)
    wn = 100 .+ randn(rng, 60, 60)
    wmask = rand(rng, 60, 60) .< 0.3
    wmask[[1:7; 54:60], :] .= false
    wmask[:, [1:7; 54:60]] .= false
    μw, dw = proc_continuous(wn, wmask; Np = 7, widx = 21, ndraw = 2, ftype = 64)
    @test 0.9 < std([dw[:, :, 1][wmask] .- μw[wmask]; dw[:, :, 2][wmask] .- μw[wmask]]) < 1.1
end
