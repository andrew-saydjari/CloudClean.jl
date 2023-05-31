@testset "preprocess.jl" begin
    out = kstar_circle_mask(3;rlim=1)
    ref = [
        true   false   true   ;
        false  false   false  ;
        true   false   true   ;
        ]
    @test out == ref

    cx = (0:17) * ones(18)'
    cy = ones(18) * (0:17)'
    out = im_subrng(2,4,cx,cy,16,16,2,2,4,4,1,1,4,4)
    @test out[1] == 4:9
    @test out[2] == 12:17
    @test length(out[3]) == 20
    out = im_subrng(1,1,cx,cy,16,16,2,2,4,4,1,1,4,4)
    @test out[1] == 0:5
    out = im_subrng(4,4,cx,cy,15,15,2,2,4,4,1,1,4,4)
    @test out[1] == 12:16

    refin = [
        11.0  12.0   11.0 ;
        32.0  60.0   19.0 ;
        22.0  25.0   20.0 ;
        ]
    add_noise!(refin,2;seed=2021)
    ref = [
        12.0   9.5  10.0 ;
        31.5  70.5  18.0 ;
        19.0  25.5  19.5 ;
        ]
    @test abs(refin[2,2]-70.5) < 10

    ttt_testim = ones(33,33)
    ttt_bmaskim = zeros(Bool,33,33)
    ttt_bmaskim[16,16] = true
    ttt_bimage = ones(33,33)
    ttt_bimageI = ones(Int,33,33)
    ttt_testim2 = zeros(33,33)
    ttt_bmaskim2 = zeros(Bool,33,33)
    ttt_goodpix = ones(Bool,33,33)
    prelim_infill!(ttt_testim,ttt_bmaskim,ttt_bimage,ttt_bimageI,ttt_testim2,ttt_bmaskim2,ttt_goodpix,widx=19,widy=19)
    @test ttt_testim2[16,16] == 1.0
    prelim_infill!(ttt_testim,ttt_bmaskim,ttt_bimage,ttt_bimageI,ttt_testim2,ttt_bmaskim2,ttt_goodpix,widx=1,widy=1)
    @test ttt_testim2[16,16] == 1.0
    prelim_infill!(ttt_testim,ttt_bmaskim,ttt_bimage,ttt_bimageI,ttt_testim2,ttt_bmaskim2,ttt_goodpix,widx=19,widy=19,ftype=64)
    @test ttt_testim2[16,16] == 1.0

    function test_sig_iqr()
        # Test case 1
        x1 = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        expected1 = iqr(x1) / 1.34896
        @test sig_iqr(x1) ≈ expected1
    
        # Test case 2
        x2 = [3, 7, 2, 1, 5, 8, 4, 9, 6]
        expected2 = iqr(x2) / 1.34896
        @test sig_iqr(x2) ≈ expected2
    
        # Add more test cases if needed
    
        # Test case with an empty array
        x_empty = []
        @test_throws DomainError sig_iqr(x_empty)
    end
    test_sig_iqr()

    function test_add_sky_noise()
        # Test case
        testim2_0 = [1.0, 2.0, 3.0, 4.0, 5.0]
        testim2_1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        maskim = [true, false, true, true, false]
        sig_iqr = 0.5
        add_sky_noise!(testim2_1, maskim, sig_iqr)
    
        @test testim2_1[.!maskim] == testim2_0[.!maskim]
        @test maximum(abs.(testim2_1.-testim2_0)).<2
    
        # Test case with an empty array
        testim2_empty = []
        maskim_empty = []
        sig_iqr_empty = 0.1
        @test_throws DomainError add_sky_noise!(testim2_empty, maskim_empty, sig_iqr_empty)
    end
    test_add_sky_noise()



end