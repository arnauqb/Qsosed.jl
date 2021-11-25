using Qsosed, Test


@testset "Test corona properties" begin
    # comparison with Xspec model

    bh = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    corona =
        Corona(bh, hard_xray_fraction = 0.02, electron_energy = 100, reprocessing = true)
    Ledd = compute_eddington_luminosity(bh)

    computed_radius = compute_corona_radius(corona)
    @test computed_radius ≈ 8.9008961989970832 rtol = 1e-3

    temp_corona = compute_disk_temperature(bh, computed_radius)
    @test temp_corona ≈ 107102.98423666063 rtol = 1e-1

    dissip_lumin = compute_corona_dissipated_luminosity(corona) / Ledd
    @test dissip_lumin ≈ 2.0183853632627917e-2 rtol = 1e-2

    corona_luminosity = compute_corona_luminosity(corona) / Ledd
    @test corona_luminosity ≈ 3.6390444796974425e-2 rtol = 1e-1

    corona_photon_index = compute_corona_photon_index(corona)
    @test corona_photon_index ≈ 2.2826826043284854 rtol = 1e-2

    # Second black hole

    bh = BlackHole(1e6 * M_SUN, 0.025, 0.0)
    @test gravity_radius(bh) ≈ 1161.0642639605840 rtol = 1e-2
    corona =
        Corona(bh, hard_xray_fraction = 0.02, electron_energy = 100, reprocessing = true)
    Ledd = compute_eddington_luminosity(bh)

    computed_radius = compute_corona_radius(corona)
    @test computed_radius ≈ 89.731710842601927  rtol = 1e-2

    temp_corona = compute_disk_temperature(bh, computed_radius)
    @test temp_corona ≈ 47643.787634112188 rtol = 1e-1

    dissip_lumin = compute_corona_dissipated_luminosity(corona) / Ledd
    @test dissip_lumin ≈ 2.0001856311471888e-2 rtol = 1e-2

    corona_luminosity = compute_corona_luminosity(corona) / Ledd
    @test corona_luminosity ≈ 2.0621184497821025e-2 rtol = 1e-1

    corona_photon_index = compute_corona_photon_index(corona)
    @test corona_photon_index ≈ 1.6483970267857697  rtol = 1e-1

end
