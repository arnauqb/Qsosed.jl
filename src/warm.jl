export Warm, compute_warm_photon_index, compute_ywarm, compute_warm_spectral_luminosity
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

function compute_warm_spectral_luminosity(corona, warm, energy_range)
    gamma = warm.photon_index
    kt_e = warm.electron_energy
    kt_w = compute_disk_temperature(corona.bh, warm.radius) * K_B
    kt_w_kev = kt_w * ERG_TO_KEV
    params = [gamma, kt_e, kt_w_kev, 0, 0]
    photons_per_bin = compute_compton_photons_per_bin(energy_range, params)
    mid_energies = (energy_range[1:(end - 1)] + energy_range[2:end]) / 2
    lumin_per_bin = vcat([0], photons_per_bin[2:end] .* mid_energies)
    total_lumin = sum(lumin_per_bin)
    target_luminosity = compute_disk_luminosity(corona.bh, corona.radius, warm.radius)
    ret = lumin_per_bin * (target_luminosity / total_lumin)
    ret = ret[2:end] ./ diff(energy_range)
    return vcat([0], ret)
end
