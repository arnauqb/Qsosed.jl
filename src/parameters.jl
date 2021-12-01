export Parameters
using YAML

struct Parameters{T<:AbstractFloat,S<:Flag,U<:Bool,V<:Int,W<:String}
    # T float, S flag, U bool, V int, W string
    M::T
    mdot::T
    spin::T
    mu_nucleon::T
    reprocessing::U
    hard_xray_fraction::T
    corona_electron_energy::T
    warm_electron_energy::T
    warm_photon_index::T
    reflection_albedo::T
    min_energy::T
    max_energy::T
    uv_min_energy::T
    uv_max_energy::T
    xray_min_energy::T
    xray_max_energy::T
    function Parameters(;
        M=1e8,
        mdot=0.5,
        spin=0.0,
        mu_nucleon=0.61,
        reprocessing=true,
        hard_xray_fraction=0.02,
        corona_electron_energy=100.0,
        warm_electron_energy=0.2,
        warm_photon_index=2.5,
        reflection_albedo=0.3,
        min_energy=1e-4, # kev
        max_energy=200.0,
        uv_min_energy=0.00387,
        uv_max_energy=0.06,
        xray_min_energy=0.06,
        xray_max_energy=1.0,
    )
        return new{Float64, Flag, Bool, Int, String}(
            M,
            mdot,
            spin,
            mu_nucleon,
            reprocessing,
            hard_xray_fraction,
            corona_electron_energy,
            warm_electron_energy,
            warm_photon_index,
            reflection_albedo,
            min_energy,
            max_energy,
            uv_min_energy,
            uv_max_energy,
            xray_min_energy,
            xray_max_energy,
        )
    end
end


function Parameters(file_path::String)
    config = YAML.load_file(file_path, dicttype=Dict{Symbol, Any})
    return Parameters(;config...)
end
