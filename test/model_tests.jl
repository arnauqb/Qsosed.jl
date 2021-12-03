using Qsosed, Test

@testset "read parameters" begin
    config_path = String(@__DIR__) * "/../configs/config_example.yaml"
    model = QsosedModel(config_path)
    @test model.corona.bh == model.bh
    @test model.bh.M == 1e8 * M_SUN
    @test model.bh.mdot == 0.5
    @test model.bh.spin == 0.0
    @test model.corona.reprocessing == true
    @test model.corona.hard_xray_fraction == 0.02
    @test model.corona.electron_energy == 100
    @test model.warm.electron_energy == 0.2
    @test model.warm.photon_index == 2.5
end
