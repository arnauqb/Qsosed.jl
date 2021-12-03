using Qsosed, Test

@testset "Compton test" begin
    parameters = [2.5, 100, 50, 0, 0]
    ear = 10 .^ range(-4, 2, length=50)
    @test compute_compton_photons_per_bin(ear, parameters) == donthcomp(ear, parameters)
end
