@testset "per_loc.jl" begin

    function test_gen_pix_mask_trivial()
        # Test case 1
        kmasked2d = falses(33, 33)
        kstar, kcond = gen_pix_mask_trivial(kmasked2d, Np=33)
        @test kstar == kmasked2d
        @test kcond == 0
    
        # Test case 2
        kmasked2d = trues(33, 33)
        kstar, kcond = gen_pix_mask_trivial(kmasked2d, Np=33)
        @test kstar == kmasked2d
        @test kcond == 0
    
        # Test case 3
        kmasked2d = rand(Bool, 33, 33)
        kstar, kcond = gen_pix_mask_trivial(kmasked2d, Np=33)
        @test kstar == kmasked2d
        @test kcond == Np^2 - count(kstar)
    end
    test_gen_pix_mask_trivial()


    function test_gen_pix_mask_circ()
        # Test case 1
        kmasked2d = falses(33, 33)
        circmask = ones(Bool, 33, 33)
        kstar, kcond = gen_pix_mask_circ(kmasked2d, circmask, Np=33)
        @test kstar == kmasked2d
        @test kcond == 0
    
        # Test case 2
        kmasked2d = trues(33, 33)
        circmask = falses(33, 33)
        kstar, kcond = gen_pix_mask_circ(kmasked2d, circmask, Np=33)
        @test kstar == kmasked2d
        @test kcond == 33^2
    
        # Test case 3
        kmasked2d = rand(Bool, 33, 33)
        circmask = trues(33, 33)
        kstar, kcond = gen_pix_mask_circ(kmasked2d, circmask, Np=33)
        @test kstar == (kmasked2d .| circmask)
        @test kcond == 33^2 - count(kstar)
    end
    test_gen_pix_mask_circ()

    function test_build_cov()
        # Test case 1
        cov = zeros(Float64, 100, 100)
        μ = zeros(Float64, 100)
        cx = 50
        cy = 50
        bimage = ones(Float64, 100, 100)
        bism = ones(Float64, 100, 100, 100, 100)
        Np = 9
        widx = 3
        widy = 3
        build_cov!(cov, μ, cx, cy, bimage, bism, Np, widx, widy)
    
        # Assert some conditions on the computed covariance matrix
        @test size(cov) == (100, 100)
        @test sum(cov) == 0
        @test cov[50, 50] ≈ 1.0
    end
    test_build_cov()

    function test_condCovEst_wdiag()
        # Test case 1
        cov_loc = rand(100, 100)
        μ = rand(100)
        km = rand(Bool, 100)
        data_in = rand(100)
        Np = 33
        export_mean = true
        n_draw = 10
        seed = 2022
    
        out = condCovEst_wdiag(cov_loc, μ, km, data_in, Np=Np, export_mean=export_mean, n_draw=n_draw, seed=seed)
    
        # Assert conditions on the output
        @test length(out) == 2
        @test size(out[1]) == (100,)
        @test size(out[2]) == (100, n_draw)

    end
    test_condCovEst_wdiag()
end