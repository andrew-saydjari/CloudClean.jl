using Random
using Distributions

export prelim_infill!
export kstar_circle_mask
export im_subrng
export add_sky_noise! # Document
export sig_iqr # Document

function sig_iqr(x)
    return iqr(x)/1.34896
end

function add_sky_noise!(testim2,maskim,sig_iqr;seed=2021)
    rng = MersenneTwister(seed)
    dist = Distributions.Normal(0,sig_iqr)
    for i in eachindex(testim2)
        if maskim[i]
            testim2[i] += rand(rng, dist)
        end
    end
end

"""
    kstar_circle_mask(Np;rlim=256) -> circmask

    Generates a Bool mask for pixels beyond a given (squared) radius of the center of an image.

    # Arguments:
    - `Np`: size of image stamp

    # Keywords:
    - `rlim`: squared radius (in pixels^2) beyond which pixels should be masked (default 256)

    # Outputs:
    - `circmask`: static Bool mask used for assigning pixels beyond some radius of the stellar center as "ignored"
"""
function kstar_circle_mask(Np;rlim=256)
    halfNp = (Np-1) ÷ 2
    x = (-halfNp:halfNp)' .* ones(Int,Np)
    y = ones(Int,Np)' .* (-halfNp:halfNp)
    R = x.^2 .+ y.^2
    return R.>rlim
end

"""
    im_subrng(jx,jy,cx,cy,sx,sy,px0,py0,stepx,stepy,padx,pady,tilex,tiley) -> xrng, yrng, star_ind

    Computes the flux a star must have so that the PSF-based masking using `thr`
    would require a larger stamp area. Used for computational savings.

    # Arguments:
    - `jx`: tile index along x
    - `jy`: tile index along y
    - `cx`: list of stellar x-coordinates
    - `cy`: list of stellar y-coordinates
    - `sx`: size of original image in x
    - `sy`: size of original image in y
    - `px0`: maximal padding in x to account for stars outside image
    - `py0`: maximal padding in y to account for stars outside image
    - `stepx`: tiling step size in x
    - `stepy`: tiling step size in y
    - `padx`: tile padding in x required to account for local stamp size, sample size, and pixels outside the image
    - `pady`: tile padding in x required to account for local stamp size, sample size, and pixels outside the image
    - `tilex`: total number of tile divisions along x
    - `tiley`: total number of tile divisions along y

    # Outputs:
    - `xrng`: slicing range of the tile in x
    - `yrng`: slicing range of the tile in y
    - `star_ind`: Bool mask of all stars falling within the tile (subimage)
"""
function im_subrng(jx,jy,cx,cy,sx,sy,px0,py0,stepx,stepy,padx,pady,tilex,tiley)
    lowbx = (1 + (jx-1)*stepx)
    uppbx = (1 + jx*stepx-1)
    lowby = (1 + (jy-1)*stepy)
    uppby = (1 + jy*stepy-1)
    # pad right and top by 2 to fix ragged
    if uppbx > sx
        uppbx = sx
    end
    if uppby > sy
        uppby = sy
    end
    xrng = (lowbx-padx):(uppbx+padx)
    yrng = (lowby-pady):(uppby+pady)

    if jx == 1
        lowbx-=px0
    end
    if jx == tilex
        uppbx+=px0
    end
    if jy == 1
        lowby-=py0
    end
    if jy == tiley
        uppby+=py0
    end

    star_ind = findall((lowbx-0.5 .< cx .<= uppbx+0.5) .& (lowby-0.5 .< cy .<= uppby+0.5))

    return xrng, yrng, star_ind
end

"""
    prelim_infill!(testim,bmaskim,bimage,bimageI,testim2,bmaskim2,goodpix,ccd;widx=19,widy=19,ftype::Int=32,widmult=1.4)

    This intial infill replaces masked pixels with a guess based on a smoothed
    boxcar. For large masked regions, the smoothing scale is increased. If this
    iteration takes too long/requires too strong of masking, the masked pixels
    are replaced with the median of the image.

    We use 3 copies of the input image and mask image. The first
    is a view (with reflective boundary condition padding) with the pixels to be infilled
    replaced with zeros, the second is allocated to hold various smoothings of the image,
    and the third holds the output image which contains our best infill guess. A final Bool
    array of size corresponding to the image is used to keep track of pixels that have safe
    infill values.

    # Arguments:
    - `testim`: input image which requires infilling
    - `bmaskim`: input mask indicating which pixels require infilling
    - `bimage`: preallocated array for smoothed version of input image
    - `bimageI`: preallocated array for smoothed mask counting the samples for each estimate
    - `testim2`: inplace modified ouptut array for infilled version of image input
    - `bmaskim2`: inplace modified mask to keep track of which pixels still need infilling
    - `goodpix`: preallocated array for Bool indexing pixels with good infill
    - `ccd`: string name of FITS extension for verbose cmdline printing

    # Keywords:
    - `widx`: initial size of boxcar smoothing window in x (default 19)
    - `widy`: initial size of boxcar smoothing window in y (default 19)
    - `ftype::Int`: determine the Float precision, 32 is Float32, otherwise Float64
    - `widmult`: multiplicative factor for increasing the smoothing scale at each iteration step
"""
function prelim_infill!(testim,bmaskim,bimage,bimageI,testim2,bmaskim2,goodpix,ccd;widx=19,widy=19,ftype::Int=32,widmult=1.4)
    if ftype == 32
        T = Float32
    else
        T = Float64
    end
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    (sx, sy) = size(testim)

    widxMax = round(Int,(widmult^10)*(widx-1)/2)*2+1
    widyMax = round(Int,(widmult^10)*(widy-1)/2)*2+1
    ΔxMax = (widxMax-1)÷2
    ΔyMax = (widyMax-1)÷2

    #the masked entries in testim must be set to 0 so they drop out of the mean
    testim[bmaskim] .= 0;
    bmaskim2 .= copy(bmaskim)
    testim2 .= copy(testim)

    #hopefully replace with the reflected indexedx arrays
    in_image = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(ΔxMax,ΔxMax)));
    in_mask = ImageFiltering.padarray(.!bmaskim,ImageFiltering.Pad(:reflect,(ΔyMax,ΔyMax)));

    #loop to try masking at larger and larger smoothing to infill large holes
    cnt=0
    while any(bmaskim2) .& (cnt .< 10)
        # this double deep view should be only 1 deep ideally... need the internal unwrap
        @views in_image1 = in_image[(1-Δx):(sx+Δx),(1-Δy):(sy+Δy)]
        @views in_mask1 = in_mask[(1-Δx):(sx+Δx),(1-Δy):(sy+Δy)]
        (sx1, sy1) = size(in_image1)
        tot = zeros(T,sx1)
        totI = zeros(Int,sx1)
        boxsmooth!(bimage,in_image1,tot,widx,widy)
        boxsmooth!(bimageI,in_mask1,totI,widx,widy)

        goodpix .= (bimageI .> 10)

        testim2[bmaskim2 .& goodpix] .= (bimage./bimageI)[bmaskim2 .& goodpix]
        bmaskim2[goodpix] .= false

        # update loop params
        cnt+=1
        widx*=1.4
        widy*=1.4
        widx = round(Int,(widx-1)/2)*2+1
        widy = round(Int,(widy-1)/2)*2+1
        Δx = (widx-1)÷2
        Δy = (widy-1)÷2
    end
    println("Infilling $ccd completed after $cnt rounds with final width (widx,widy) = ($widx,$widy)")
    flush(stdout)

    #catastrophic failure fallback
    if cnt == 10
        testim2[bmaskim2] .= StatsBase.median(testim)
        println("Infilling Failed Badly")
        flush(stdout)
    end
    return
end