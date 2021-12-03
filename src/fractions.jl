export total_xray_fraction,
    radial_uv_fraction,
    disk_uv_fraction,
    spectral_band_fraction_frequency,
    spectral_band_fraction

"""
    spectral_band_fraction_frequency(bb, low, high)

Computes the amount of black body radiance that is emitted in
a particular frequency band.
"""
function spectral_band_fraction_frequency(bb::BlackBody, low, high)
    total_radiance = radiance(bb)
    band_radiance = spectral_band_radiance_frequency(bb, low, high)
    return band_radiance / total_radiance
end
"""
    spectral_band_fraction(bb, low, high)
Computes the amount of black body radiance that is emitted in
a particular energy band.
"""

function spectral_band_fraction(bb::BlackBody, low, high)
    total_radiance = radiance(bb)
    band_radiance = spectral_band_radiance(bb, low, high)
    return band_radiance / total_radiance
end

"""
    disk_uv_fraction(bh::BlackHole, r; uv_low_kev = 0.00387, uv_high_kev = 0.06)
Computes the fraction of luminosity in the UV band, for a ring at disk r.
The UV band is defined as uv_low_kev -- uv_high_kev.
"""
function disk_uv_fraction(bh::BlackHole, r; uv_low_kev = 0.00387, uv_high_kev = 0.06)
    if r <= bh.isco
        return 0.0
    end
    temperature = compute_disk_temperature(bh, r)
    bb = BlackBody(temperature)
    return spectral_band_fraction(bb, uv_low_kev, uv_high_kev)
end

"""
    warm_uv_fraction(warm, bh, r; uv_low_kev = 0.00387, uv_high_kev = 0.06)
Computes the fraction of luminosity in the UV band, for a ring at disk r.
The UV band is defined as uv_low_kev -- uv_high_kev.
"""
function compute_warm_uv_fraction(
    corona,
    warm,
    bh;
    uv_low_kev = 0.00387,
    uv_high_kev = 0.06,
    energy_low = 1e-5,
    energy_high = 1e3,
    n_energy = 5000,
)
    energy_range = 10 .^ range(log10(energy_low), log10(energy_high), length = n_energy)
    warm_spectral_luminosity = compute_warm_spectral_luminosity(corona, warm, energy_range)
    uv_lower_b = searchsortedfirst(energy_range, uv_low_kev)
    uv_upper_b = searchsortedfirst(energy_range, uv_high_kev)
    delta_energy_range = diff(energy_range)
    uv_lumin = sum(
        warm_spectral_luminosity[uv_lower_b:uv_upper_b] .*
        delta_energy_range[uv_lower_b:uv_upper_b],
    )
    total_lumin = sum(warm_spectral_luminosity[2:end] .* delta_energy_range)
    return uv_lumin / total_lumin
end

"""
    disk_xray_fraction(bh::BlackHole, r; xray_low_kev = 0.00387, xray_high_kev = 0.06)
Computes the fraction of luminosity in the UV band, for a ring at disk r.
The UV band is defined as uv_low_kev -- uv_high_kev.
"""
function disk_xray_fraction(bh::BlackHole; xray_low_kev = 0.06, xray_high_kev = 1)
    total_lumin = compute_bolometric_luminosity(bh)
    xray_lumin = disk_spectral_band_radiance(bh, xray_low_kev, xray_high_kev)
    xray_lumin / total_lumin
end


"""
    total_xray_fraction(corona, warm; xray_low_kev = 0.00387, xray_high_kev = 0.06)
Computes the fraction of luminosity that is in the X-ray fraction.
"""
function total_xray_fraction(
    corona,
    warm;
    energy_low = 1e-5,
    energy_high = 1e3,
    xray_low_kev = 0.06,
    xray_high_kev = 1,
    n_energy = 5000,
)
    energy_range = 10 .^ range(log10(energy_low), log10(energy_high), length = n_energy)
    cl = compute_corona_spectral_luminosity(corona, warm, energy_range)
    dl = Qsosed.compute_disk_spectral_luminosity(
        corona.bh,
        energy_range,
        r_min = warm.radius,
        r_max = gravity_radius(corona.bh),
    )
    wl = compute_warm_spectral_luminosity(corona, warm, energy_range)
    total_spectral_lumin = cl + dl + wl
    xray_lower_b = searchsortedfirst(energy_range, xray_low_kev)
    xray_upper_b = searchsortedfirst(energy_range, xray_high_kev)
    delta_energy_range = diff(energy_range)
    xray_lumin = sum(
        total_spectral_lumin[xray_lower_b:xray_upper_b] .*
        delta_energy_range[xray_lower_b:xray_upper_b],
    )
    total_lumin = sum(total_spectral_lumin[2:end] .* delta_energy_range)
    return xray_lumin / total_lumin
end
total_xray_fraction(
    model;
    energy_low = 1e-5,
    energy_high = 1e3,
    xray_low_kev = 0.06,
    xray_high_kev = 1,
    n_energy = 5000,
) = total_xray_fraction(
    model.corona,
    model.warm,
    energy_low = energy_low,
    energy_high = energy_high,
    xray_low_kev = xray_low_kev,
    xray_high_kev = xray_high_kev,
    n_energy = n_energy,
)

"""
    radial_uv_fraction(corona, warm; uv_low_kev = 0.00387, uv_high_kev = 0.06)
Fraction of power emitted in the UV band per radius.
"""
function radial_uv_fraction(
    corona,
    warm;
    uv_low_kev = 0.00387,
    uv_high_kev = 0.06,
    n_r = 5000,
)
    r_range =
        10 .^ range(log10(corona.bh.isco), log10(gravity_radius(corona.bh)), length = n_r)
    warm_uv_fraction = compute_warm_uv_fraction(
        corona,
        warm,
        corona.bh,
        uv_low_kev = uv_low_kev,
        uv_high_kev = uv_high_kev,
    )
    ret = Float64[]
    for r in r_range
        if r <= corona.radius
            push!(ret, 0)
        elseif corona.radius < r < warm.radius
            push!(ret, warm_uv_fraction)
        else
            push!(
                ret,
                disk_uv_fraction(
                    corona.bh,
                    r,
                    uv_low_kev = uv_low_kev,
                    uv_high_kev = uv_high_kev,
                ),
            )
        end
    end
    return r_range, ret
end
radial_uv_fraction(model; n_r = 5000) = radial_uv_fraction(
    model.corona,
    model.warm,
    uv_low_kev = model.parameters.uv_min_energy,
    uv_high_kev = model.parameters.uv_max_energy,
    n_r = n_r,
)

