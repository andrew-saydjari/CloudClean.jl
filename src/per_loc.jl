using LinearAlgebra

export build_cov!
export condCovEst_wdiag
export gen_pix_mask_trivial
export gen_pix_mask_circ

"""
    gen_pix_mask_trivial(kmasked2d; Np=33) -> kstar, kcond

Flatten a pixel mask and calculate the number of pixels used for the conditional infill.

# Arguments:
- `kmasked2d`: A 2D array representing the masked pixels.

# Keywords:
- `Np`: An optional integer specifying the number of pixels in a side (default: 33).

# Returns:
- `kstar`: A flattened version of the input `kmasked2d` array.
- `kcond`: The count of unmasked pixels in the `kstar` array.

# Examples:
```julia
julia> kmasked2d = rand(Bool, 33, 33)
julia> kstar, kcond = gen_pix_mask_trivial(kmasked2d, Np=33)
```
"""
function gen_pix_mask_trivial(kmasked2d;Np=33)
    kstar = kmasked2d
    kcond = Np^2-count(kstar)

    return kstar[:], kcond
end

"""  
    gen_pix_mask_circ(kmasked2d, circmask; Np=33) -> kstar, kcond

Generate a circular pixel mask and calculate the number of pixels used for the conditional infill.

# Arguments:
- `kmasked2d`: A 2D array representing the masked pixels.
- `circmask`: A 2D array representing the circular mask.

# Keywords:
- `Np`: An optional integer specifying the number of pixels in a side (default: 33).

# Returns:
- `kstar`: A copy of the input `kmasked2d` array with circular masking applied.
- `kcond`: The count of unmasked pixels in the `kstar` array.

# Examples:
```julia
julia> kmasked2d = rand(Bool, 33, 33)
julia> circmask = kstar_circle_mask(33,rlim=256)
julia> kstar, kcond = gen_pix_mask_circ(kmasked2d, circmask, Np=33)
    ```
"""
function gen_pix_mask_circ(kmasked2d,circmask;Np=33)
    kstar = kmasked2d .| circmask
    kcond = Np^2-count(kstar)

    return kstar[:], kcond
end

"""
    condCovEst_wdiag(cov_loc,μ,km,data_in;Np=33,export_mean=false,n_draw=0) -> out

Using a local covariance matrix estimate `cov_loc` and a set of known ("good") pixels `km`, this function computes a prediction for the mean value of masked pixels and the covariance matrix of the masked pixels. The output list can conditionally include the mean reconstruction and draws from the distribution of reconstructions.

# Arguments:
- `cov_loc`: local covariance matrix
- `μ`: vector containing mean value for each pixel in the patch
- `km`: unmasked pixels
- `data_in`: input image

# Keywords:
- `Np`: size of local covariance matrix in pixels (default 33)
- `export_mean`: when true, returns the mean conditional prediction for the "hidden" pixels (default false)
- `n_draw`: when nonzero, returns that number of realizations of the conditional infilling (default 0)

# Outputs:
- `out[1]`: input image returned with masked pixels replaced with mean prediction
- `out[2]`: input image returned with masked pixels replaced with a draw from the predicted distribution
"""
function condCovEst_wdiag(cov_loc,μ,km,data_in;Np=33,export_mean=false,n_draw=0,seed=2022)
    k = .!km
    kstar = km
    cov_kk = Symmetric(cov_loc[k,k])
    cov_kkstar = cov_loc[k,kstar];
    cov_kstarkstar = cov_loc[kstar,kstar];
    icov_kkC = cholesky(cov_kk)
    icovkkCcovkkstar = icov_kkC\cov_kkstar
    predcovar = Symmetric(cov_kstarkstar - (cov_kkstar'*icovkkCcovkkstar))
    ipcovC = cholesky(predcovar)

    @views uncond_input = data_in[:]
    @views cond_input = data_in[:].- μ

    kstarpredn = (cond_input[k]'*icovkkCcovkkstar)'
    kstarpred = kstarpredn .+ μ[kstar]

    out = []
    if export_mean
        mean_out = copy(data_in)
        mean_out[kstar] .= kstarpred
        push!(out,mean_out)
    end
    if n_draw != 0
        sqrt_cov = ipcovC.U
        noise = randn(n_draw,size(sqrt_cov)[1])*sqrt_cov

        draw_out = repeat(copy(data_in)[:],outer=[1 n_draw])
        draw_out[kstar,:] .= repeat(kstarpred,outer=[1 n_draw]) .+ noise'
        push!(out,draw_out)
    end

    return out
end

"""
    build_cov!(cov::Array{T,2},μ::Array{T,1},cx::Int,cy::Int,bimage::Array{T,2},bism::Array{T,4},Np::Int,widx::Int,widy::Int) where T <:Union{Float32,Float64}

Constructs the local covariance matrix and mean for an image patch of size `Np` x `Np` pixels around a location
of interest (`cx`,`cy`). The construction is just a lookup of pixel values from the stored boxcar-smoothed copies
of the input image times itself shifted in `bism`. Passing the smoothed image `bimage` and the widths of the boxcar
mean `widx` and `widy` is helpful for the mean and normalization. The covariance and mean are updated in place
for speed since this operation may be performed billions of times since we construct a new covariance matrix for
every detection. Math may either be performed `Float32` or `Float64`.

# Arguments:
- `cov::Array{T,2}`: preallocated output array for local covariance matrix
- `μ::Array{T,1}`: preallocated output vector for local mean
- `cx::Int`: x-coordinate of the center of the local region
- `cy::Int`: y-coordinate of the center of the local region
- `bimage::Array{T,2}`: boxcar smoothed unshifted image
- `bism::Array{T,4}`: boxcar-smoothed image products for all shifts
- `Np::Int`: size of local covariance matrix in pixels
- `widx::Int`: width of boxcar window in x which determines size of region used for samples for the local covariance estimate
- `widy::Int`: width of boxcar window in y which determines size of region used for samples for the local covariance estimate
"""
function build_cov!(cov::Array{T,2},μ::Array{T,1},cx::Int,cy::Int,bimage::Array{T,2},bism::Array{T,4},Np::Int,widx::Int,widy::Int) where T <:Union{Float32,Float64}
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    halfNp = (Np-1) ÷ 2
    Δr, Δc = cx-(halfNp+1), cy-(halfNp+1)
    # Δr, Δc = cx-(halfNp-1), cy-(halfNp-1)
    for dc=0:Np-1       # column shift loop
        pcr = 1:Np-dc
        for dr=1-Np:Np-1# row loop, incl negatives
            if (dr < 0) & (dc == 0)
                continue
            end
            if dr >= 0
                prr = 1:Np-dr
            end
            if (dr < 0) & (dc > 0)
                prr = 1-dr:Np
            end

            for pc=pcr, pr=prr
                i = ((pc   -1)*Np)+pr
                j = ((pc+dc-1)*Np)+pr+dr
                @inbounds μ1μ2 = bimage[pr+Δr,pc+Δc]*bimage[pr+dr+Δr,pc+dc+Δc]/((widx*widy)^2)
                @inbounds t = bism[pr+Δr,pc+Δc,dr+Np,dc+1]/(widx*widy) - μ1μ2
                @inbounds cov[i,j] = t
                @inbounds cov[j,i] = t
                if i == j
                    @inbounds μ[i] = μ1μ2
                end
            end
        end
    end
    cov .*= (widx*widy)/((widx*widy)-1)
    return
end