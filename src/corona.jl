export Corona,
    compute_corona_covering_factor,
    compute_corona_seed_luminosity,
    compute_corona_dissipated_luminosity,
    compute_corona_photon_index,
    compute_corona_radius,
    compute_corona_luminosity

using QuadGK, Roots

struct Corona{T<:Float64,U<:Bool}
    bh::BlackHole
    hard_xray_fraction::T
    electron_energy::T
    radius::T
    reprocessing::U
end
function Corona(bh; hard_xray_fraction, electron_energy, reprocessing)
    radius = compute_corona_radius(bh, hard_xray_fraction)
    return Corona(bh, Float64(hard_xray_fraction), Float64(electron_energy), Float64(radius), reprocessing)
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
compute_corona_dissipated_luminosity(corona) = compute_corona_dissipated_luminosity(corona.bh, corona.hard_xray_fraction)

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
        diff = truncated_disk_lumin - compute_corona_dissipated_luminosity(bh, hard_xray_fraction)
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
    return compute_corona_seed_luminosity(corona) + compute_corona_dissipated_luminosity(corona)
end
