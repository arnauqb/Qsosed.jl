using QuadGK
export BlackBody,
    spectral_radiance,
    spectral_radiance_frequency,
    spectral_band_radiance,
    spectral_band_radiance_frequency,
    radiance,
    disk_nt_rel_factors,
    gravity_radius,
    disk_flux,
    compute_disk_temperature,
    compute_disk_luminosity,
    compute_disk_spectral_luminosity,
    compute_disk_spectral_flux,
    compute_disk_spectral_flux_at_radius,
    compute_disk_flux_at_radius

struct BlackBody
    T::Float64
end

kev_to_hz(energy) = energy / ERG_TO_KEV / H_PLANCK
lambda_peak(T) = 0.2897771955 / T
nu_peak(T) = C / lambda_peak(T)
nu_peak(bb::BlackBody) = nu_peak(bb.T)

"""
Spectral radiance of a black body bb for a frequency. Input frequency is
in Hz, and output units are erg / (s cm^2 sr Hz)
"""
function spectral_radiance_frequency(bb::BlackBody, frequency)
    a = 2 * H_PLANCK * frequency^3 / C^2
    b = 1 / (exp((H_PLANCK * frequency) / (K_B * bb.T)) - 1)
    return a * b
end

"""
Spectral radiance of a black body at a given energy. The input energy
has units of keV, and the result output is in units of
keV / (s cm² keV sr)
"""
function spectral_radiance(bb::BlackBody, energy)
    frequency = kev_to_hz(energy)
    rad = spectral_radiance_frequency(bb, frequency)
    return rad
end

"""
    spectral_band_radiance_frequency(bb::BlackBody, low, high)
Computes the radiance on a frequency band from low to high. low and high
have to be in Hz, return units are erg / sr / m^2
"""
function spectral_band_radiance_frequency(bb::BlackBody, low, high)
    low = max(low, nu_peak(bb) / 1e4)
    high = min(high, nu_peak(bb) * 1e4)
    integral, err = quadgk(x -> spectral_radiance_frequency(bb, x), low, high, rtol = 1e-8)
    return integral
end

"""
Computes the radiance on an energetic band from low to high.
"""
function spectral_band_radiance(bb::BlackBody, low, high)
    low = kev_to_hz(low)
    high = kev_to_hz(high)
    low = max(low, nu_peak(bb) / 1e3)
    high = min(high, nu_peak(bb) * 1e3)
    integral, err = quadgk(x -> spectral_radiance_frequency(bb, x), low, high, rtol = 1e-8)
    return integral
end

radiance(bb::BlackBody) = SIGMA_SB * bb.T^4 / π

"""
Novikov-Thorne relativistic factors for the AD spectrum.
If the passed radius is smaller than ISCO, it returns 0.

Parameters
----------
radius
    disc radius measured in R_g
spin
    normalized black hole spin (-1, 1)
isco
    innermost stable circular orbit measured in Rg
"""
function disk_nt_rel_factors(radius, spin, isco)
    if radius <= isco
        return 0.0
    end
    yms = sqrt(isco)
    y1 = 2 * cos((acos(spin) - pi) / 3)
    y2 = 2 * cos((acos(spin) + pi) / 3)
    y3 = -2 * cos(acos(spin) / 3)
    y = sqrt(radius)
    C = 1 - 3 / radius + 2 * spin / radius^1.5
    B = 3 * (y1 - spin)^2 * log((y - y1) / (yms - y1)) / (y * y1 * (y1 - y2) * (y1 - y3))
    B += 3 * (y2 - spin)^2 * log((y - y2) / (yms - y2)) / (y * y2 * (y2 - y1) * (y2 - y3))
    B += 3 * (y3 - spin)^2 * log((y - y3) / (yms - y3)) / (y * y3 * (y3 - y1) * (y3 - y2))
    A = 1 - yms / y - 3 * spin * log(y / yms) / (2 * y)
    factor = (A - B) / C
    return factor
end

disk_nt_rel_factors(bh, r) = disk_nt_rel_factors(r, bh.spin, bh.isco)

"""
Self-gravity radius as described by Laor & Netzer (1989).
"""
function gravity_radius(bh)
    mass = (bh.M / (1e9 * M_SUN))
    alpha = 0.1 # assumption
    r_sg = 2150 * mass^(-2 / 9) * bh.mdot^(4 / 9) * alpha^(2 / 9)# / bh.Rg
    return r_sg
end


"""
Computes the disk flux of the accretion disc of black hole ``bh``
at a given radius ``r`` (in Rg). See equation 5 of https://arxiv.org/pdf/2001.04720.pdf
"""
function disk_flux(bh, r)
    Mdot = compute_mass_accretion_rate(bh)
    NT_factors = disk_nt_rel_factors(bh, r)
    return 3 * G * bh.M * Mdot / (8 * π * (r * bh.Rg)^3) * NT_factors
end

"""
Computes the temperature of the disk at a radius r (in units of Rg).
"""
function compute_disk_temperature(bh, r)
    flux = disk_flux(bh, r)
    return (flux / SIGMA_SB)^(1 / 4)
end

"""
    compute_disk_luminosity(bh, r_min, r_max)
computes the luminosity of the disk between r_min and r_max
"""
function compute_disk_luminosity(bh, r_min, r_max)
    constant = SIGMA_SB * 4 * π * bh.Rg^2
    integ, _ = quadgk(
        r -> r * compute_disk_temperature(bh, r)^4,
        r_min,
        r_max,
        atol = 0,
        rtol = 1e-2,
    )
    return integ * constant
end

"""
    compute_disk_spectral_flux_at_radius(bh, radius, energy_range)
Returns an array with the photon flux at each bin. (# photons / cm^2)
"""
function compute_disk_spectral_flux_at_radius(bh, radius, energy)
    ret = Float64[]
    bb = BlackBody(compute_disk_temperature(bh, radius))
    freq = kev_to_hz(energy)
    radiance = spectral_radiance_frequency(bb, freq)
    return radiance * π / HZ_TO_KEV
end

"""
    compute_disk_flux_at_radius(bh, radius, energy_range)
Returns an array with the photon flux at each bin. (# photons / cm^2)
"""
function compute_disk_flux_at_radius(bh, radius, low = 1e-5, high = 1e2)
    ret = Float64[]
    bb = BlackBody(compute_disk_temperature(bh, radius))
    radiance = spectral_band_radiance(bb, low, high)
    return radiance * π
end

function compute_disk_spectral_luminosity(bh, energy_range; r_min=6, r_max=1600)
    ret = Float64[]
    for energy in energy_range
        integ, _ = quadgk(
            r -> r * compute_disk_spectral_flux_at_radius(bh, r, energy),
            r_min,
            r_max,
            rtol = 1e-8,
            atol = 0,
        )
        push!(ret, integ)
    end
    return ret * 2 * 2π * bh.Rg^2
end
