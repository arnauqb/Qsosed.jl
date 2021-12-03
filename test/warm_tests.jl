using Qsosed, Test

@testset "Test warm component" begin

    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    corona =
        Corona(bh, hard_xray_fraction = 0.02, electron_energy = 100, reprocessing = true)
    warm = Warm(corona, electron_energy = 0.2, photon_index = 2.5)
    @test compute_warm_photon_index(warm) == 2.5
    @test warm.electron_energy == 0.2
    @test warm.radius ≈ 17.801792397994166 rtol = 1e-2

    temp_warm = compute_disk_temperature(bh, warm.radius)
    @test temp_warm ≈ 87745.726703792374 rtol = 1e-1
end
