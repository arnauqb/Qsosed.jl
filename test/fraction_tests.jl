using Test, Qsosed

@testset "Disk UV fraction" begin
    bb = BlackBody(1e4)
    @test spectral_band_fraction_frequency(bb, 1e12, 1e18) ≈ 1 rtol = 1e-8
    @test spectral_band_fraction(bb, 1e-5, 1e2) ≈ 1 rtol = 1e-6

    bh = BlackHole(1e8 * M_SUN, 0.5, 0)
    @test disk_uv_fraction(bh, 1.0) ≈ 0
    uvf1 = disk_uv_fraction(bh, 6.01)
    uvf2 = disk_uv_fraction(bh, 30)
    uvf3 = disk_uv_fraction(bh, 200)
    @test uvf1 < uvf2
    @test uvf3 < uvf2
end

@testset "Radial UV fractions" begin
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    corona =
        Corona(bh, hard_xray_fraction = 0.02, electron_energy = 100, reprocessing = true)
    warm = Warm(corona, electron_energy = 0.2, photon_index = 2.5)
    rr, uvf = radial_uv_fraction(corona, warm, n_r = 500)
    @test length(rr) == length(uvf) == 500
    for (i, r) in enumerate(rr)
        if r <= corona.radius
            @test uvf[i] == 0.0
            continue
        end
        if r <= warm.radius
            @test uvf[i] < disk_uv_fraction(bh, r)
            continue
        end
        if r > warm.radius
            @test uvf[i] == disk_uv_fraction(bh, r)
        end
    end
end

@testset "X-Ray fraction" begin
    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    corona =
        Corona(bh, hard_xray_fraction = 0.02, electron_energy = 100, reprocessing = true)
    warm = Warm(corona, electron_energy = 0.2, photon_index = 2.5)
    xf = total_xray_fraction(corona, warm)
    @test 0.1 < xf < 0.2
end
