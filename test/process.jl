@testset "process.jl" begin

    function test_proc_continuous()
        # Test case 1
        raw_image = rand(100, 100)
        mask_image = rand(Bool, 100, 100)
        Np = 33
        widx = 129
        widy = widx
        tilex = 1
        tiley = tilex
        seed = 2021
        ftype = 32
        ndraw = 0
    
        out = proc_continuous(raw_image, mask_image, Np=Np, widx=widx, widy=widy, tilex=tilex, tiley=tiley, seed=seed, ftype=ftype, ndraw=ndraw)
    
        # Assert conditions on the output
        @test size(out) == (100, 100)
    end
    test_proc_continuous()

    function test_proc_discrete()
        # Test case 1
        x_locs = [50]
        y_locs = [50]
        raw_image = rand(100, 100)
        mask_image = rand(Bool, 100, 100)
        Np = 33
        widx = 129
        widy = widx
        tilex = 1
        tiley = tilex
        seed = 2021
        ftype = 32
        ndraw = 0
    
        out = proc_discrete(raw_image, mask_image, Np=Np, widx=widx, widy=widy, tilex=tilex, tiley=tiley, seed=seed, ftype=ftype, ndraw=ndraw)
    
        # Assert conditions on the output
        @test size(out) == (100, 100)
    end
    test_proc_discrete()
end