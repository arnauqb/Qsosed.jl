export compute_compton_photon_flux 

function compute_compton_photon_flux(energy_range, parameters)
    photon_number_per_bin = donthcomp(energy_range, parameters)
    photon_number_flux = photon_number_per_bin ./ energy_range
    return photon_number_flux
end
