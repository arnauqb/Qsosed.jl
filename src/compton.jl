export compute_compton_photon_flux, compute_compton_photons_per_bin

function compute_compton_photons_per_bin(energy_range, parameters)
    photon_number_per_bin = donthcomp(energy_range, parameters)
    return photon_number_per_bin
end

