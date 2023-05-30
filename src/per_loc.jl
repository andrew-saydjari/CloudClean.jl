using LinearAlgebra

export build_cov!
export condCovEst_wdiag
export condCovEst_wdiag_full
export gen_pix_mask_AKS

function gen_pix_mask_trivial(kmasked2d;Np=33)
    kstar = kmasked2d
    kcond = Np^2-count(kstar)

    return kstar[:], kcond
end

"""
    condCovEst_wdiag(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=33,export_mean=false,n_draw=0,diag_on=true) -> out

    Using a local covariance matrix estimate `cov_loc` and a set of known ("good") pixels `km`
    and "hidden" pixels `kpsf2d`, this function computes a prediction for the mean value
    of the `kpsf2d` pixels and the covariance matrix of the `kpsf2d` pixels. In terms of
    statistics use to adjust the photometry of a star, we are only interested in the
    pixels masked as a result of the star (i.e. not a detector defect or cosmic ray nearby).
    The residual image `data_in` and a model of the counts above the background coming from the
    star `stars_in` for the local patch are also inputs of the function. Correction factors for
    the photometric flux and flux uncertainities are outputs as well as a chi2 value for the
    "good" pixels. The output list can conditionally include the mean reconstruction and
    draws from the distribution of reconstructions.

    # Arguments:
    - `cov_loc`: local covariance matrix
    - `μ`: vector containing mean value for each pixel in the patch
    - `km`: unmasked pixels
    - `kpsf2d`: pixels masked due to the star of interest
    - `data_in`: (non-infilled) residual image in local patch
    - `psft`: static array (image) of the stellar PSF

    # Keywords:
    - `Np`: size of local covariance matrix in pixels (default 33)
    - `export_mean`: when true, returns the mean conditional prediction for the "hidden" pixels (default false)
    - `n_draw`: when nonzero, returns that number of realizations of the conditional infilling (default 0)
    - `diag_on`: flag for adding to the pixelwise uncertainty based on the photoelectron counts of the modeled star (default true)

    # Outputs:
    - `out[1][1]`: flux uncertainity of the star
    - `out[1][2]`: flux uncertainity of the star assuming the covariance matrix were diagonal
    - `out[1][3]`: flux correction which must be added to correct the input flux estimate
    - `out[1][4]`: flux correction coming from the residuals (fdb_res)
    - `out[1][5]`: flux correction coming from the predicted background (fdb_pred)
    - `out[1][6]`: chi2 for the "good" pixels under `cov_loc` as a metric on how good our assumptions are
    - `out[2]`: local region (image) with "hidden" pixels replaced by the mean conditional estimate (optional output)
    - `out[end]`: local region (image) with "hidden" pixels replaced by the draws from the conditional distribution (optional output). Array is flattened to npix x n_draw.
"""
function condCovEst_wdiag_discrete(cov_loc,μ,km,kpsf2d,data_in,stars_in,psft;Np=33,export_mean=false,n_draw=0,diag_on=true)
    k = .!km
    kstar = kpsf2d[:]
    if diag_on
        for i=1:Np*Np cov_loc[i,i] += stars_in[i] end
    end
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
    @views p = psft[kpsf2d][:]
    ipcovCp = ipcovC\p

    #@views std_wdiag = sqrt(abs(sum((pw.^(2)).*diag(predcovar[kpsf1d_kstar,kpsf1d_kstar]))))/sum(p2w)
    @views var_wdb = (p'*ipcovCp)

    var_diag = 0
    for i in 1:size(predcovar)[1]
        var_diag += (p[i]^2)/predcovar[i,i]
    end

    @views resid_mean = (uncond_input[kstar]'*ipcovCp)./var_wdb
    @views pred_mean = (kstarpred'*ipcovCp)./var_wdb

    #if we can afford it, a nice check would be how good of a covariance matrix
    #cov is for the k pixels (which we think are clean)
    chi2 = cond_input[k]'*(icov_kkC\cond_input[k])

    # Currently limited to the Np region. Often useful to have some context with a larger
    # surrounding region... TO DO to implement
    out = []
    push!(out,[sqrt(var_wdb^(-1)) sqrt(var_diag^(-1)) pred_mean-resid_mean resid_mean pred_mean chi2])
    if export_mean
        mean_out = copy(data_in)
        mean_out[kstar] .= kstarpred
        push!(out,mean_out)
    end
    if n_draw != 0
        covsvd = svd(predcovar)
        sqrt_cov = covsvd.V*diagm(sqrt.(covsvd.S))*covsvd.Vt;
        noise = sqrt_cov*randn(size(sqrt_cov)[1],n_draw)

        draw_out = repeat(copy(data_in)[:],outer=[1 n_draw])
        draw_out[kstar,:] .= repeat(kstarpred,outer=[1 n_draw]) .+ noise
        push!(out,draw_out)
    end

    return out
end

function condCovEst_wdiag_continuous(cov_loc,μ,km,data_in;Np=33,export_mean=false,n_draw=0,seed=2022)
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

    # Currently limited to the Np region. Often useful to have some context with a larger
    # surrounding region... TO DO to implement
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