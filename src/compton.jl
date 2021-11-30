export compute_compton_photon_flux, compute_compton_photons_per_bin

function compute_compton_photons_per_bin(energy_range, parameters)
    photon_number_per_bin = donthcomp(energy_range, parameters)
    return photon_number_per_bin
end

function compute_compton_photon_flux(energy_range, parameters)
    photon_number_per_bin = donthcomp(energy_range, parameters)
    photon_number_flux = photon_number_per_bin[2:end] ./ diff(energy_range)
    return vcat([0.0], photon_number_flux)
end
