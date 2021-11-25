export Warm, compute_warm_photon_index
struct Warm{T<:Float64}
    electron_energy::T
    photon_index::T
    radius::T
end

function Warm(corona::Corona; electron_energy, photon_index)
    radius = 2 * corona.radius
    return Warm(electron_energy, photon_index, radius)
end

compute_warm_photon_index(warm) = warm.photon_index
