using Qsosed, Test

@testset "Test that they agree with ADT implementation in python" begin
    ear = 10 .^ range(-4, 2, length=50)
    #parameters = [2.27963222215913, 100, 0.01676681194689287, 0, 0]
    for gamma in range(2, 10, length=10)
        for kt_e in range(10, 1000, length=10)
            for t_corona_kev in range(0.001, 0.1, length=10)
                parameters = [gamma, kt_e, t_corona_kev, 0, 0]
                expected = donthcomp_py(ear, parameters)
                results = donthcomp(ear, parameters)
                @test expected ≈ results rtol = 1e-3
            end
        end
    end
end
