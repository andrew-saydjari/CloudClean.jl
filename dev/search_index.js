var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"api/#Covariance-Construction-Functions","page":"API Reference","title":"Covariance Construction Functions","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"cov_avg!\nboxsmooth!\noutest_bounds","category":"page"},{"location":"api/#CloudClean.cov_avg!","page":"API Reference","title":"CloudClean.cov_avg!","text":"cov_avg!(bimage, ism, bism, in_image; Np::Int=33, widx::Int=129, widy::Int=129, ftype::Int=32)\n\nKey function for constructing the (shifted and multiplied) versions of the input image used to quickly\nestimate the local covariance matrix at a large number of locations. The main output is in the preallocated\n`bism` which is used as an input to `build_cov!`.\n\n# Arguments:\n- `bimage`: preallocated output array for the boxcar smoothed unshifted image\n- `ism`: preallocated intermediate array for the input image times itself shifted\n- `bism`: preallocated output array to store boxcar-smoothed image products for all shifts\n- `in_image`: input image the local covariance of which we want to estimate\n\n# Keywords:\n- `Np::Int`: size of local covariance matrix in pixels (default 33)\n- `widx::Int`: width of boxcar window in x which determines size of region used for samples for the local covariance estimate (default 129)\n- `widy::Int`: width of boxcar window in y which determines size of region used for samples for the local covariance estimate (default 129)\n- `ftype::Int`: determine the Float precision, 32 is Float32, otherwise Float64\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.boxsmooth!","page":"API Reference","title":"CloudClean.boxsmooth!","text":"boxsmooth!(out::AbstractArray, arr::AbstractArray, tot::Array{T,1}, widx::Int, widy::Int)\n\nBoxcar smooths an input image (or paddedview) `arr` with window size `widx` by\n`widy`. We pass the original image size `sx` and `sy` to help handle image views.\n\n# Arguments:\n- `out::AbstractArray`: preallocated output array for the boxcar smoothed image\n- `arr::AbstractArray`: input array for which boxcar smoothing is computed (generally paddedview)\n- `tot::Array{T,1}`: preallocated array to hold moving sums along 1 dimension\n- `widx::Int`: size of boxcar smoothing window in x\n- `widy::Int`: size of boxcar smoothing window in y\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.outest_bounds","page":"API Reference","title":"CloudClean.outest_bounds","text":"outest_bounds(cx,sx) -> px0\n\nHelper function to find maximum padding in pixels required to accomodate all query points `cx` outside of the image size 1:`sx`.\n\n# Arguments:\n- `cx`: list of integer star centers (in either x or y)\n- `sx`: image dimension along the axis indexed by `cx`\n\n# Outputs:\n- `px0`: maximum padding in pixels required to accomodate all query points\n\n\n\n\n\n","category":"function"},{"location":"api/#Per-Location-Functions","page":"API Reference","title":"Per Location Functions","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"gen_pix_mask_trivial\ngen_pix_mask_circ\ncondCovEst_wdiag\nbuild_cov!","category":"page"},{"location":"api/#CloudClean.gen_pix_mask_trivial","page":"API Reference","title":"CloudClean.gen_pix_mask_trivial","text":"    gen_pix_mask_trivial(kmasked2d; Np=33) -> kstar, kcond\n\nFlatten a pixel mask and calculate the number of pixels used for the conditional infill.\n\n# Arguments:\n- `kmasked2d`: A 2D array representing the masked pixels.\n\n# Keywords:\n- `Np`: An optional integer specifying the number of pixels in a side (default: 33).\n\n# Returns:\n- `kstar`: A flattened version of the input `kmasked2d` array.\n- `kcond`: The count of unmasked pixels in the `kstar` array.\n\n# Examples:\n```julia\njulia> kmasked2d = rand(Bool, 33, 33)\njulia> kstar, kcond = gen_pix_mask_trivial(kmasked2d, Np=33)\n```\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.gen_pix_mask_circ","page":"API Reference","title":"CloudClean.gen_pix_mask_circ","text":"    gen_pix_mask_circ(kmasked2d, circmask; Np=33) -> kstar, kcond\n\n<!-- -->\nGenerate a circular pixel mask and calculate the number of pixels used for the conditional infill.\n\n# Arguments:\n- `kmasked2d`: A 2D array representing the masked pixels.\n- `circmask`: A 2D array representing the circular mask.\n\n# Keywords:\n- `Np`: An optional integer specifying the number of pixels in a side (default: 33).\n\n# Returns:\n- `kstar`: A copy of the input `kmasked2d` array with circular masking applied.\n- `kcond`: The count of unmasked pixels in the `kstar` array.\n\n# Examples:\n```julia\njulia> kmasked2d = rand(Bool, 33, 33)\njulia> circmask = kstar_circle_mask(33,rlim=256)\njulia> kstar, kcond = gen_pix_mask_circ(kmasked2d, circmask, Np=33)\n```\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.condCovEst_wdiag","page":"API Reference","title":"CloudClean.condCovEst_wdiag","text":"    condCovEst_wdiag(cov_loc,μ,km,data_in;Np=33,export_mean=false,n_draw=0) -> out\n\nUsing a local covariance matrix estimate `cov_loc` and a set of known (\"good\") pixels `km`, this function computes a prediction for the mean value of masked pixels and the covariance matrix of the masked pixels. The output list can conditionally include the mean reconstruction and draws from the distribution of reconstructions.\n\n# Arguments:\n- `cov_loc`: local covariance matrix\n- `μ`: vector containing mean value for each pixel in the patch\n- `km`: unmasked pixels\n- `data_in`: input image\n\n# Keywords:\n- `Np`: size of local covariance matrix in pixels (default 33)\n- `export_mean`: when true, returns the mean conditional prediction for the \"hidden\" pixels (default false)\n- `n_draw`: when nonzero, returns that number of realizations of the conditional infilling (default 0)\n\n# Outputs:\n- `out[1]`: input image returned with masked pixels replaced with mean prediction\n- `out[2]`: input image returned with masked pixels replaced with a draw from the predicted distribution\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.build_cov!","page":"API Reference","title":"CloudClean.build_cov!","text":"    build_cov!(cov::Array{T,2},μ::Array{T,1},cx::Int,cy::Int,bimage::Array{T,2},bism::Array{T,4},Np::Int,widx::Int,widy::Int) where T <:Union{Float32,Float64}\n\nConstructs the local covariance matrix and mean for an image patch of size `Np` x `Np` pixels around a location\nof interest (`cx`,`cy`). The construction is just a lookup of pixel values from the stored boxcar-smoothed copies\nof the input image times itself shifted in `bism`. Passing the smoothed image `bimage` and the widths of the boxcar\nmean `widx` and `widy` is helpful for the mean and normalization. The covariance and mean are updated in place\nfor speed since this operation may be performed billions of times since we construct a new covariance matrix for\nevery detection. Math may either be performed `Float32` or `Float64`.\n\n# Arguments:\n- `cov::Array{T,2}`: preallocated output array for local covariance matrix\n- `μ::Array{T,1}`: preallocated output vector for local mean\n- `cx::Int`: x-coordinate of the center of the local region\n- `cy::Int`: y-coordinate of the center of the local region\n- `bimage::Array{T,2}`: boxcar smoothed unshifted image\n- `bism::Array{T,4}`: boxcar-smoothed image products for all shifts\n- `Np::Int`: size of local covariance matrix in pixels\n- `widx::Int`: width of boxcar window in x which determines size of region used for samples for the local covariance estimate\n- `widy::Int`: width of boxcar window in y which determines size of region used for samples for the local covariance estimate\n\n\n\n\n\n","category":"function"},{"location":"api/#Image-Preprocessing","page":"API Reference","title":"Image Preprocessing","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"prelim_infill!\nsig_iqr\nim_subrng\nadd_sky_noise!\nadd_noise!\nkstar_circle_mask","category":"page"},{"location":"api/#CloudClean.prelim_infill!","page":"API Reference","title":"CloudClean.prelim_infill!","text":"    prelim_infill!(testim,bmaskim,bimage,bimageI,testim2,bmaskim2,goodpix;widx=19,widy=19,ftype::Int=32,widmult=1.4)\n\nThis intial infill replaces masked pixels with a guess based on a smoothed\nboxcar. For large masked regions, the smoothing scale is increased. If this\niteration takes too long/requires too strong of masking, the masked pixels\nare replaced with the median of the image.\n\nWe use 3 copies of the input image and mask image. The first\nis a view (with reflective boundary condition padding) with the pixels to be infilled\nreplaced with zeros, the second is allocated to hold various smoothings of the image,\nand the third holds the output image which contains our best infill guess. A final Bool\narray of size corresponding to the image is used to keep track of pixels that have safe\ninfill values.\n\n# Arguments:\n- `testim`: input image which requires infilling\n- `bmaskim`: input mask indicating which pixels require infilling\n- `bimage`: preallocated array for smoothed version of input image\n- `bimageI`: preallocated array for smoothed mask counting the samples for each estimate\n- `testim2`: inplace modified ouptut array for infilled version of image input\n- `bmaskim2`: inplace modified mask to keep track of which pixels still need infilling\n- `goodpix`: preallocated array for Bool indexing pixels with good infill\n\n# Keywords:\n- `widx`: initial size of boxcar smoothing window in x (default 19)\n- `widy`: initial size of boxcar smoothing window in y (default 19)\n- `ftype::Int`: determine the Float precision, 32 is Float32, otherwise Float64\n- `widmult`: multiplicative factor for increasing the smoothing scale at each iteration step\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.sig_iqr","page":"API Reference","title":"CloudClean.sig_iqr","text":"    sig_iqr(x)\n\nCalculate the normalized interquartile range (IQR) of an array as a robust measure of standard deviation.\n\n# Arguments\n- `x`: A 1D array or iterable.\n\n# Returns\n- The normalized IQR, computed as the IQR divided by 1.34896.\n\n# Examples\n```julia\njulia> x = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]\njulia> result = sig_iqr(x)\n```\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.im_subrng","page":"API Reference","title":"CloudClean.im_subrng","text":"    im_subrng(jx,jy,cx,cy,sx,sy,px0,py0,stepx,stepy,padx,pady,tilex,tiley) -> xrng, yrng, star_ind\n\nComputes the flux a star must have so that the PSF-based masking using `thr`\nwould require a larger stamp area. Used for computational savings.\n\n# Arguments:\n- `jx`: tile index along x\n- `jy`: tile index along y\n- `cx`: list of stellar x-coordinates\n- `cy`: list of stellar y-coordinates\n- `sx`: size of original image in x\n- `sy`: size of original image in y\n- `px0`: maximal padding in x to account for stars outside image\n- `py0`: maximal padding in y to account for stars outside image\n- `stepx`: tiling step size in x\n- `stepy`: tiling step size in y\n- `padx`: tile padding in x required to account for local stamp size, sample size, and pixels outside the image\n- `pady`: tile padding in x required to account for local stamp size, sample size, and pixels outside the image\n- `tilex`: total number of tile divisions along x\n- `tiley`: total number of tile divisions along y\n\n# Outputs:\n- `xrng`: slicing range of the tile in x\n- `yrng`: slicing range of the tile in y\n- `star_ind`: Bool mask of all stars falling within the tile (subimage)\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.add_sky_noise!","page":"API Reference","title":"CloudClean.add_sky_noise!","text":"    add_sky_noise!(testim2, maskim, sig_iqr; seed=2021)\n\nAdd sky noise to pixels in an image specified by a given mask.\n\n# Arguments:\n- `testim2`: A mutable 2D array representing the image.\n- `maskim`: A 2D array representing the mask.\n- `sig_iqr`: The standard deviation of the noise distribution, generally calculated as the normalized IQR.\n\n# Keywords:\n- `seed`: An optional integer specifying the random number generator seed (default: 2021).\n\n# Returns:\n- Modifies `testim2` in place by adding sky noise to the masked pixels.\n\n# Examples:\n```julia\njulia> testim2 = rand(100, 100)\njulia> maskim = rand(Bool, 100, 100)\njulia> sig_iqr = 0.5\njulia> add_sky_noise!(testim2, maskim, sig_iqr, seed=2021)\n```\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.add_noise!","page":"API Reference","title":"CloudClean.add_noise!","text":"    add_noise!(testim2,gain;seed=2021)\n\nAdds noise to an image that matches the Poisson noise of the pixel counts.\nA random seed to set a local random generator is provided for reproducible unit testing.\n\n# Arguments:\n- `testim2`: input image which had infilling\n- `gain`: gain of detector to convert from photon count noise to detector noise\n\n# Keywords:\n- `seed`: random seed for random generator\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.kstar_circle_mask","page":"API Reference","title":"CloudClean.kstar_circle_mask","text":"    kstar_circle_mask(Np;rlim=256) -> circmask\n\nGenerates a Bool mask for pixels beyond a given (squared) radius of the center of an image.\n\n# Arguments:\n- `Np`: size of image stamp\n\n# Keywords:\n- `rlim`: squared radius (in pixels^2) beyond which pixels should be masked (default 256)\n\n# Outputs:\n- `circmask`: static Bool mask used for assigning pixels beyond some radius of the stellar center as \"ignored\"\n\n\n\n\n\n","category":"function"},{"location":"api/#Full-Image-Processing-Functions","page":"API Reference","title":"Full Image Processing Functions","text":"","category":"section"},{"location":"api/","page":"API Reference","title":"API Reference","text":"proc_discrete\nproc_continuous","category":"page"},{"location":"api/#CloudClean.proc_discrete","page":"API Reference","title":"CloudClean.proc_discrete","text":"proc_discrete(x_locs, y_locs, raw_image, mask_image; Np=33, widx=129, widy=widx, tilex=1, tiley=tilex, seed=2021, ftype::Int=32, ndraw=0) -> out_mean, out_draw\n\nProcess an image with a mask, replacing masked pixels with either a mean or draw from a distribution resembling the local pixel-pixel covariance structure in the image.\n\nArguments:\n\nx_locs: A 1D array representing the location centers (in the x coordinate) for infilling.\ny_locs: A 1D array representing the location centers (in the y coordinate) for infilling.\nraw_image: A 2D array representing the input image.\nmask_image: A 2D array representing the mask.\n\nKeywords:\n\nNp: An optional integer specifying the number of pixels in a side (default: 33).\nwidx: An optional integer specifying the width of the region used for training the local covariance in the x-direction (default: 129).\nwidy: An optional integer specifying the width of the region used for training the local covariance in the y-direction (default: widx).\ntilex: An optional integer specifying the number of tiles in the x-direction for subdividing the image (default: 1).\ntiley: An optional integer specifying the number of tiles in the y-direction for subdividing the image (default: tilex).\nseed: An optional integer specifying the random number generator seed (default: 2021).\nftype: An optional integer specifying the floating-point precision type (32 or 64) (default: 32).\nrlim: Radius limit for the radial mask beyond which pixels are not used for conditioning (units are pixels^2). (default: 625)\nndraw: An optional integer specifying the number of draws of samples from the statistical distribution of possible masked pixel values (default: 0).\n\nReturns\n\nIf ndraw is 0, returns the debiased image as a 2D array.\nIf ndraw is greater than 0, returns the debiased image as a 2D array and an array of ndraw draws.\n\nExamples\n\njulia> raw_image = rand(100, 100)\njulia> mask_image = kstar_circle_mask(100,rlim=256)\njulia> result = proc_continuous([50],[50],raw_image, mask_image, Np=33, widx=129, seed=2021)\n\n\n\n\n\n","category":"function"},{"location":"api/#CloudClean.proc_continuous","page":"API Reference","title":"CloudClean.proc_continuous","text":"    proc_continuous(raw_image, mask_image; Np=33, widx=129, widy=widx, tilex=1, tiley=tilex, seed=2021, ftype::Int=32, ndraw=0) -> out_mean, out_draw\n\nProcess an image with a mask, replacing masked pixels with either a mean or draw from a distribution resembling the local pixel-pixel covariance structure in the image.\n\n# Arguments:\n- `raw_image`: A 2D array representing the input image.\n- `mask_image`: A 2D array representing the mask.\n\n# Keywords:\n- `Np`: An optional integer specifying the number of pixels in a side (default: 33).\n- `widx`: An optional integer specifying the width of the region used for training the local covariance in the x-direction (default: 129).\n- `widy`: An optional integer specifying the width of the region used for training the local covariance in the y-direction (default: widx).\n- `tilex`: An optional integer specifying the number of tiles in the x-direction for subdividing the image (default: 1).\n- `tiley`: An optional integer specifying the number of tiles in the y-direction for subdividing the image (default: tilex).\n- `seed`: An optional integer specifying the random number generator seed (default: 2021).\n- `ftype`: An optional integer specifying the floating-point precision type (32 or 64) (default: 32).\n- `ndraw`: An optional integer specifying the number of draws of samples from the statistical distribution of possible masked pixel values (default: 0).\n\n# Returns\n- If `ndraw` is 0, returns the debiased image as a 2D array.\n- If `ndraw` is greater than 0, returns the debiased image as a 2D array and an array of `ndraw` draws.\n\n# Examples\n```julia\njulia> raw_image = rand(100, 100)\njulia> mask_image = rand(Bool, 100, 100)\njulia> result = proc_continuous(raw_image, mask_image, Np=33, widx=129, seed=2021)\n```\n\n\n\n\n\n","category":"function"},{"location":"#CloudClean.jl","page":"Introduction","title":"CloudClean.jl","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"(Image: GitHub) (Image: Build Status) (Image: Coverage Status)","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"Image infilling algorithm with focus on statistically accuracy","category":"page"},{"location":"#Installation","page":"Introduction","title":"Installation","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"CloudClean is a registered package so a stable version can be installed using Pkg.add.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.add(\"CloudClean\")","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"For the most recent development version, install directly from the GitHub","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"import Pkg\nPkg.add(url=\"https://github.com/andrew-saydjari/CloudClean.jl\")","category":"page"},{"location":"#Description","page":"Introduction","title":"Description","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"By leveraging the local pixel-pixel covarariance structure in an image, CloudClean attempts to predict the values of missing pixels based on near-by unmasked pixels. The user needs only to provide an image and a mask of \"bad\" pixels and the choice of a few hyper parameters. CloudClean has two main operating modes:","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"proc_discrete, which attempts to fill in masked data in a subregion centered on a discrete list of input points\nproc_continuous, which infills arbitrarily shaped and centered masks. The latter is like \"Photoshop\" for images with correlated structure.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"This code is based heavily upon (and is in some sense a simplification of) CloudCovErr.jl.","category":"page"},{"location":"#Example","page":"Introduction","title":"Example","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"An example of the quality of the prediction for missing pixels is demonstrated on this image from the WISE 12 μm dust map. More examples, notebooks, and documentation are in process.","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"[![\"WISE infill example\"][infill-img]][infill-url]","category":"page"},{"location":"#Contributing-and-Questions","page":"Introduction","title":"Contributing and Questions","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"This is a new piece of software. Filing an issue to report a bug or request a feature is extremely valuable in helping us prioritize what to work on, so don't hesitate.","category":"page"},{"location":"#Table-of-Contents","page":"Introduction","title":"Table of Contents","text":"","category":"section"},{"location":"","page":"Introduction","title":"Introduction","text":"Pages = [\"index.md\",\"api.md\"]","category":"page"},{"location":"","page":"Introduction","title":"Introduction","text":"<!– URLS –> [infill-img]: docs/src/assets/infillradiuswhite.gif [infill-url]: https://faun.rc.fas.harvard.edu/saydjari/CloudCovErr/thr_test.mp4","category":"page"}]
}
