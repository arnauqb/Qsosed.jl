using Qsosed
using Test

@testset "Test Black Hole properties" begin
    black_hole = BlackHole(1e8 * M_SUN, 0.5, 0.5)
    @test compute_eddington_luminosity(black_hole) ≈ 1.25750e46 * 1.17 rtol = 1e-3 #erg/s
    @test compute_bolometric_luminosity(black_hole) ≈
          0.5 * compute_eddington_luminosity(black_hole) rtol = 1e-3
    #@test mass_accretion_rate(black_hole) ≈ 1.39916e26 rtol = 1e-4
    @test black_hole.Rg ≈ 14766250380501.246 rtol = 1e-4
    @test black_hole.isco ≈ 4.2330025 rtol = 1e-4
    black_hole2 = BlackHole(1e8 * M_SUN, 0.5, 0.0)
    @test black_hole2.isco ≈ 6.0 rtol = 1e-4
    @test black_hole2.efficiency ≈ 0.057190958417936644
    @test compute_mass_accretion_rate(black_hole2) ≈ 1.2228101097715833e+26 * 1.17
end

@testset "Initialise black holes" begin
    parameters = Parameters(M=1e9, mdot=0.2, spin=0.6)
    bh = BlackHole(parameters)
    @test bh.M == 1e9 * M_SUN
    @test bh.mdot == 0.2
    @test bh.spin == 0.6
end
