using Test, Qsosed



@testset "test parse parameters" begin
    config_path = String(@__DIR__) * "/../configs/config_example.yaml"
    parameters = Parameters(config_path)
    @test parameters.M == 1e8
    @test parameters.mdot == 0.5
    @test parameters.mdot == 0.5
    @test parameters.spin == 0.0
    @test parameters.mu_nucleon == 0.61
    @test parameters.reprocessing == true
    @test parameters.hard_xray_fraction == 0.02
    @test parameters.corona_electron_energy == 100
    @test parameters.warm_electron_energy == 0.2
    @test parameters.warm_photon_index == 2.5
    @test parameters.reflection_albedo == 0.3
    @test parameters.min_energy == 1e-4 # kev
    @test parameters.max_energy == 200
    @test parameters.uv_min_energy == 0.00387
    @test parameters.uv_max_energy == 0.06
    @test parameters.xray_min_energy == 0.06
    @test parameters.xray_max_energy == 1
end
