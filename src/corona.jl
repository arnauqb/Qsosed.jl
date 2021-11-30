export Corona,
    compute_corona_covering_factor,
    compute_corona_seed_luminosity,
    compute_corona_dissipated_luminosity,
    compute_corona_photon_index,
    compute_corona_radius,
    compute_corona_luminosity,
    compute_reprocessed_flux,
    compute_corona_photon_flux

using QuadGK, Roots

struct Corona{T<:Float64,U<:Bool}
    bh::BlackHole
    hard_xray_fraction::T
    electron_energy::T
    radius::T
    reprocessing::U
end
function Corona(bh::BlackHole; hard_xray_fraction, electron_energy, reprocessing)
    radius = compute_corona_radius(bh, hard_xray_fraction)
    return Corona(
        bh,
        Float64(hard_xray_fraction),
        Float64(electron_energy),
        Float64(radius),
        reprocessing,
    )
end
function Corona(bh::BlackHole, parameters::Parameters)
    return Corona(
        bh,
        hard_xray_fraction = parameters.hard_xray_fraction,
        electron_energy = parameters.corona_electron_energy,
        reprocessing = parameters.reprocessing,
    )
end

"""
    compute_covering_factor(corona)

Corona covering factor as seen from the disk at a radius r > corna radius.
"""
function compute_corona_covering_factor(corona, radius)
    if radius < corona.radius
        error("Radius can't be smaller than the corona.")
        return 0.0
    end
    theta0 = asin(corona.radius / radius)
    ret = theta0 - 0.5 * sin(2 * theta0)
    return ret
end

"""
    compute_corona_seed_luminosity(corona)

Seed photon luminosity intercepted from the warm region and the outer disk. Calculated assuming a truncated disk and spherical hot flow geometry.
"""
function compute_corona_seed_luminosity(corona)
    function f(radius)
        return radius *
               compute_disk_temperature(corona.bh, radius)^4 *
               compute_corona_covering_factor(corona, radius)
    end
    integ, _ = quadgk(f, corona.radius, gravity_radius(corona.bh), atol = 0, rtol = 1e-3)
    seed_lumin = 4 * corona.bh.Rg^2 * SIGMA_SB * integ
    return seed_lumin
end

"""
    compute_corona_dissipated_luminosity(corona)

Intrinsic luminosity from the Corona. This is assumed to be a constant fraction of the Eddington luminosity, regardless of actual accretion rate.
"""
function compute_corona_dissipated_luminosity(bh, hard_xray_fraction)
    return hard_xray_fraction * compute_eddington_luminosity(bh)
end
compute_corona_dissipated_luminosity(corona) =
    compute_corona_dissipated_luminosity(corona.bh, corona.hard_xray_fraction)

"""
    compute_corona_photon_index
Photon index (Gamma) for the corona SED. The functional form is assumed to be L_nu = k nu ^(-alpha) = k nu^( 1 - gamma ), where alpha = gamma - 1.
Computed using equation 14 of Beloborodov (1999).
"""
function compute_corona_photon_index(corona)
    gamma_cor =
        7 / 3 *
        (
            compute_corona_dissipated_luminosity(corona) /
            compute_corona_seed_luminosity(corona)
        )^(-0.1)
    return gamma_cor
end

"""
    compute_corona_radius
Finds the radius of the corona
"""
function compute_corona_radius(bh, hard_xray_fraction)
    function corona_radius_kernel(bh, radius)
        truncated_disk_lumin = compute_disk_luminosity(bh, bh.isco, radius)
        diff =
            truncated_disk_lumin -
            compute_corona_dissipated_luminosity(bh, hard_xray_fraction)
        return diff
    end
    r_cor = find_zero(
        r -> corona_radius_kernel(bh, r),
        (bh.isco, gravity_radius(bh)),
        Bisection(),
        atol = 0,
        rtol = 1e-3,
    )
    return r_cor
end
compute_corona_radius(corona) = corona.radius

"""
    compute_corona_luminosity(corona)
Total luminosity of the corona.
"""
function compute_corona_luminosity(corona)
    return compute_corona_seed_luminosity(corona) +
           compute_corona_dissipated_luminosity(corona)
end

"""
    compute_reprocessed_flux(corona, radius)
Computes the reprocessed flux from the corona at radius `radius`.
"""
function compute_reprocessed_flux(corona, radius; albedo = 0.3)
    Lhot = compute_corona_luminosity(corona)
    height = corona.radius * corona.bh.Rg
    radius = radius * corona.bh.Rg
    isco = corona.bh.isco * corona.bh.Rg
    Mdot = compute_mass_accretion_rate(corona.bh)
    aux1 = (3 * G * corona.bh.M * Mdot) / (8 * π * radius^3)
    aux2 = 2 * Lhot / (Mdot * C^2)
    aux3 = height / isco * (1 - albedo)
    aux4 = (1 + (height / radius)^2)^(-1.5)
    return aux1 * aux2 * aux3 * aux4
end


function compute_corona_photon_flux(corona::Corona, warm, energy_range, distance=1e10)
    gamma = compute_corona_photon_index(corona)
    kt_e = corona.electron_energy
    kt_c = compute_disk_temperature(corona.bh, corona.radius) * K_B
    kt_c_kev = kt_c * ERG_TO_KEV
    ywarm = compute_ywarm(warm)
    params = [gamma, kt_e, kt_c_kev * exp(ywarm), 0, 0]
    compton_flux = compute_compton_photon_flux(energy_range, params)
    println(compton_flux)
    total_compton_flux = sum(energy_range[2:end] .* diff(energy_range) .* compton_flux[2:end]) / ERG_TO_KEV
    println(total_compton_flux)
    target_flux = compute_corona_luminosity(corona) / (4 * π * distance^2)
    ret = (target_flux / total_compton_flux) * compton_flux
    return ret
end

compute_corona_photon_flux(model::Model, energy_range, distance=1e10) =
    compute_corona_photon_flux(model.corona, model.warm, energy_range, distance)

