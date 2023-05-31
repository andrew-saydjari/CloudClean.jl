using ProgressMeter
using StatsBase

export proc_continuous
export proc_discrete

function proc_continuous(raw_image,mask_image;Np=33,widx=129,widy=widx,tilex=1,tiley=tilex,seed=2021,ftype::Int=32,ndraw=0)
    radNp = (Np-1)÷2
    if ftype == 32
        T = Float32
    else
        T = Float64
    end

    # renaming to match conventions
    ref_im = raw_image
    bmaskd = mask_image
    (sx0, sy0) = size(ref_im)
    
    cartIndx = findall(bmaskd)
    x_stars = map(x->x[1],cartIndx)
    y_stars = map(x->x[2],cartIndx)

    (Nstars,) = size(x_stars)
    
    mod_im = StatsBase.median(ref_im)*ones(T,(sx0, sy0))

    testim = mod_im .- ref_im
    bimage = zeros(T,sx0,sy0)
    bimageI = zeros(Int64,sx0,sy0)
    testim2 = zeros(T,sx0,sy0)
    bmaskim2 = zeros(Bool,sx0,sy0)
    goodpix = zeros(Bool,sx0,sy0)

    prelim_infill!(testim,bmaskd,bimage,bimageI,testim2,bmaskim2,goodpix;widx=19,widy=19,ftype=ftype)
    testim .= mod_im .- ref_im #fixes current overwrite for 0 infilling

    ## calculate the star farthest outside the edge of the image in x and y
    cx = round.(Int,x_stars)
    cy = round.(Int,y_stars)
    px0 = outest_bounds(cx,sx0)
    py0 = outest_bounds(cy,sy0)

    ## these have to be allocating to get the noise model right
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    padx = Np+Δx+px0
    pady = Np+Δy+py0
    in_image = ImageFiltering.padarray(testim2,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    in_image_raw = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    in_bmaskd = ImageFiltering.padarray(bmaskd,ImageFiltering.Fill(true,(padx+2,pady+2)));
    out_mean = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    out_mean[in_bmaskd].=NaN
    out_draw = if ndraw!=1
        ImageFiltering.padarray(repeat(testim,outer=[1 1 ndraw]),ImageFiltering.Pad(:reflect,(padx+2,pady+2,0)));
    else
        ImageFiltering.padarray(repeat(testim,outer=[1 1]),ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    end
    for i=1:ndraw
        out_draw[in_bmaskd,i].=NaN
    end
    
    diffim = view(ref_im,1:sx0-1,:).-view(ref_im,2:sx0,:)
    in_sigiqr = sig_iqr(filter(.!isnan,diffim))
    
    add_sky_noise!(in_image,in_bmaskd,in_sigiqr;seed=seed)

    ## iterate over all star positions and compute errorbars/debiasing corrections
    star_stats = zeros(T,10,Nstars)
    star_k = zeros(Int32,10,Nstars)

    cov = zeros(T,Np*Np,Np*Np)
    μ = zeros(T,Np*Np)

    # some important global sizes for the loop
    cntStar0 = 0
    stepx = (sx0+2) ÷ tilex
    stepy = (sy0+2) ÷ tiley

    # precallocate the image subblocks
    in_subimage = zeros(T,stepx+2*padx,stepy+2*pady)
    ism = zeros(T,stepx+2*padx,stepy+2*pady)
    bimage = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy)
    bism = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy,2*Np-1, Np);
    for jx=1:tilex, jy=1:tiley
        xrng, yrng, star_ind = im_subrng(jx,jy,cx,cy,sx0+2,sy0+2,px0,py0,stepx,stepy,padx,pady,tilex,tiley)
        cntStar = length(star_ind)
        cntStarIter = 0
        if cntStar > 0
            in_subimage .= in_image[xrng,yrng]
            cov_avg!(bimage, ism, bism, in_subimage, widx=widx, widy=widy,Np=Np)
            offx = padx-Δx-(jx-1)*stepx
            offy = pady-Δy-(jy-1)*stepy
            for i in star_ind
                if in_bmaskd[cx[i],cy[i]]
                    #this builds in some determinism on the scan direction of the x,y
                    #processing of the image
                    build_cov!(cov,μ,cx[i]+offx,cy[i]+offy,bimage,bism,Np,widx,widy)
                    cov_stamp = cx[i]-radNp:cx[i]+radNp,cy[i]-radNp:cy[i]+radNp
                    
                    kmasked2d = in_bmaskd[cov_stamp[1],cov_stamp[2]]
                    kstar, kcond = gen_pix_mask_trivial(kmasked2d;Np=Np)
                    data_in = in_image_raw[cov_stamp[1],cov_stamp[2]]
                    try
                        stat_out = condCovEst_wdiag_continuous(cov,μ,kstar,data_in,Np=Np,export_mean=true,n_draw=ndraw,seed=seed)
                        
                        data_in[kstar].=stat_out[1][kstar]
                        in_image_raw[cov_stamp[1],cov_stamp[2]].=data_in
                        
                        data_in = out_mean[cov_stamp[1],cov_stamp[2]]
                        data_in[kstar].=stat_out[1][kstar]
                        out_mean[cov_stamp[1],cov_stamp[2]].=data_in
                        for i=1:ndraw
                            draw_in = out_draw[cov_stamp[1],cov_stamp[2],i]
                            draw_in[kstar].= stat_out[2][kstar,i]
                            out_draw[cov_stamp[1],cov_stamp[2],i].=draw_in
                        end
                        kmasked2d[kstar].=false
                        in_bmaskd[cov_stamp[1],cov_stamp[2]].=kmasked2d
                        cntStarIter += 1
                    catch
                        println("Positive Definite Error")
                    end
                end
            end
        end
        cntStar0 += cntStarIter
        println("Finished $cntStarIter of $cntStar locations in tile ($jx, $jy)")
        flush(stdout)
    end
    if ndraw>0
        return mod_im[1].-out_mean[1:sx0, 1:sy0], mod_im[1].-out_draw[1:sx0, 1:sy0, :]
    else
        return mod_im[1].-out_mean[1:sx0, 1:sy0]
    end
end

function proc_discrete(x_locs,y_locs,raw_image,mask_image;Np=33,widx=129,widy=widx,tilex=1,tiley=tilex,seed=2021,ftype::Int=32,rlim=625,ndraw=0)
    radNp = (Np-1)÷2
    if ftype == 32
        T = Float32
    else
        T = Float64
    end

    # renaming to match conventions
    ref_im = raw_image
    bmaskd = mask_image
    (sx0, sy0) = size(ref_im)

    x_stars = x_locs
    y_stars = y_locs
    (Nstars,) = size(x_stars)

    mod_im = StatsBase.median(ref_im)*ones(T,(sx0, sy0))
 
    testim = mod_im .- ref_im
    bimage = zeros(T,sx0,sy0)
    bimageI = zeros(Int64,sx0,sy0)
    testim2 = zeros(T,sx0,sy0)
    bmaskim2 = zeros(Bool,sx0,sy0)
    goodpix = zeros(Bool,sx0,sy0)

    prelim_infill!(testim,bmaskd,bimage,bimageI,testim2,bmaskim2,goodpix;widx=19,widy=19,ftype=ftype)
    testim .= mod_im .- ref_im #fixes current overwrite for 0 infilling

    ## calculate the star farthest outside the edge of the image in x and y
    cx = round.(Int,x_stars)
    cy = round.(Int,y_stars)
    px0 = outest_bounds(cx,sx0)
    py0 = outest_bounds(cy,sy0)

    ## these have to be allocating to get the noise model right
    Δx = (widx-1)÷2
    Δy = (widy-1)÷2
    padx = Np+Δx+px0
    pady = Np+Δy+py0
    in_image = ImageFiltering.padarray(testim2,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    in_image_raw = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    in_bmaskd = ImageFiltering.padarray(bmaskd,ImageFiltering.Fill(true,(padx+2,pady+2)));
    out_mean = ImageFiltering.padarray(testim,ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    out_mean[in_bmaskd].=NaN
    out_draw = if ndraw!=1
        ImageFiltering.padarray(repeat(testim,outer=[1 1 ndraw]),ImageFiltering.Pad(:reflect,(padx+2,pady+2,0)));
    else
        ImageFiltering.padarray(repeat(testim,outer=[1 1]),ImageFiltering.Pad(:reflect,(padx+2,pady+2)));
    end
    for i=1:ndraw
        out_draw[in_bmaskd,i].=NaN
    end

    diffim = view(ref_im,1:sx0-1,:).-view(ref_im,2:sx0,:)
    in_sigiqr = sig_iqr(filter(.!isnan,diffim))
    
    add_sky_noise!(in_image,in_bmaskd,in_sigiqr;seed=seed)

    ## iterate over all star positions and compute errorbars/debiasing corrections
    star_stats = zeros(T,10,Nstars)
    star_k = zeros(Int32,10,Nstars)

    cov = zeros(T,Np*Np,Np*Np)
    μ = zeros(T,Np*Np)

    # compute a radial mask for reduced num cond pixels
    circmask = kstar_circle_mask(Np,rlim=rlim)

    # some important global sizes for the loop
    cntStar0 = 0
    stepx = (sx0+2) ÷ tilex
    stepy = (sy0+2) ÷ tiley

    # precallocate the image subblocks
    in_subimage = zeros(T,stepx+2*padx,stepy+2*pady)
    ism = zeros(T,stepx+2*padx,stepy+2*pady)
    bimage = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy)
    bism = zeros(T,stepx+2*padx-2*Δx,stepy+2*pady-2*Δy,2*Np-1, Np);
    for jx=1:tilex, jy=1:tiley
        xrng, yrng, star_ind = im_subrng(jx,jy,cx,cy,sx0+2,sy0+2,px0,py0,stepx,stepy,padx,pady,tilex,tiley)
        cntStar = length(star_ind)
        if cntStar > 0
            in_subimage .= in_image[xrng,yrng]
            cov_avg!(bimage, ism, bism, in_subimage, widx=widx, widy=widy,Np=Np)
            offx = padx-Δx-(jx-1)*stepx
            offy = pady-Δy-(jy-1)*stepy
            for i in star_ind
                build_cov!(cov,μ,cx[i]+offx,cy[i]+offy,bimage,bism,Np,widx,widy)
                cov_stamp = cx[i]-radNp:cx[i]+radNp,cy[i]-radNp:cy[i]+radNp
                    
                kmasked2d = in_bmaskd[cov_stamp[1],cov_stamp[2]]
                kstar, kcond = gen_pix_mask_circ(kmasked2d,circmask;Np=Np)
                data_in = in_image_raw[cov_stamp[1],cov_stamp[2]]

                # try
                    stat_out = condCovEst_wdiag(cov,μ,kstar,data_in,Np=Np,export_mean=true,n_draw=ndraw,seed=seed)
                    
                    data_in[kstar].=stat_out[1][kstar]
                    in_image_raw[cov_stamp[1],cov_stamp[2]].=data_in
                    
                    data_in = out_mean[cov_stamp[1],cov_stamp[2]]
                    data_in[kstar].=stat_out[1][kstar]
                    out_mean[cov_stamp[1],cov_stamp[2]].=data_in
                    for i=1:ndraw
                        draw_in = out_draw[cov_stamp[1],cov_stamp[2],i]
                        draw_in[kstar].= stat_out[2][kstar,i]
                        out_draw[cov_stamp[1],cov_stamp[2],i].=draw_in
                    end
                    kmasked2d[kstar].=false
                    in_bmaskd[cov_stamp[1],cov_stamp[2]].=kmasked2d
                    cntStar0 += cntStar
                # catch
                #     println("Positive Definite Error")
                # end
            end
        end
        cntStar0 += cntStar
        println("Finished $cntStar stars in tile ($jx, $jy)")
        flush(stdout)
    end
    if ndraw>0
        return mod_im[1].-out_mean[1:sx0, 1:sy0], mod_im[1].-out_draw[1:sx0, 1:sy0, :]
    else
        return mod_im[1].-out_mean[1:sx0, 1:sy0]
    end
end