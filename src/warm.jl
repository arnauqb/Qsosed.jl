export Warm, compute_warm_photon_index, compute_ywarm
struct Warm{T<:Float64}
    electron_energy::T
    photon_index::T
    radius::T
end

function Warm(corona::Corona; electron_energy, photon_index)
    radius = 2 * corona.radius
    return Warm(electron_energy, photon_index, radius)
end
function Warm(corona, parameters)
    return Warm(
        corona,
        electron_energy = parameters.warm_electron_energy,
        photon_index = parameters.warm_photon_index,
    )
end

compute_warm_photon_index(warm) = warm.photon_index

compute_ywarm(warm) = (4 / 9 * warm.photon_index)^(-4.5)
