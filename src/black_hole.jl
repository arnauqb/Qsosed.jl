export BlackHole,
    compute_eddington_luminosity,
    compute_bolometric_luminosity,
    compute_mass_accretion_rate

struct BlackHole{T<:AbstractFloat}
    M::T
    mdot::T
    spin::T
    isco::T
    efficiency::T
    Rg::T
    function BlackHole(M, mdot, spin)
        Rg = G * M / C^2
        isco = compute_isco(spin)
        eff = compute_efficiency(spin)
        new{typeof(M)}(M, mdot, spin, isco, eff, Rg)
    end
end

function BlackHole(parameters::Parameters)
    M = parameters.M * M_SUN
    return BlackHole(M, parameters.mdot, parameters.spin)
end

"""
Computes the Eddington luminosity for the given black hole mass.
"""
function compute_eddington_luminosity(bh_mass; mu_electron = 1.17)
    return 4π * G * bh_mass * C * mu_electron / SIGMA_E
end
compute_eddington_luminosity(bh::BlackHole; mu_electron = 1.17) =
    compute_eddington_luminosity(bh.M, mu_electron = mu_electron)

"""
Computes the bolometric luminosity of the AGN given the black
hole mass and the normalized eddington rate.
"""
compute_bolometric_luminosity(bh_mass, eddington_rate; mu_electron = 1.17) =
    eddington_rate * compute_eddington_luminosity(bh_mass; mu_electron = mu_electron)
compute_bolometric_luminosity(bh::BlackHole; mu_electron = 1.17) =
    compute_bolometric_luminosity(bh.M, bh.mdot, mu_electron = 1.17)

"""
Computes the Innermost Stable Circular Orbit (ISCO) in units of Rg.

# Parameters
- spin : normalized spin parameter ∈ (-1, 1)
"""
function compute_isco(spin)
    z1 = 1 + (1 - spin^2)^(1 / 3) * ((1 + abs(spin))^(1 / 3) + (1 - abs(spin))^(1 / 3))
    z2 = sqrt(3 * spin^2 + z1^2)
    rms = 3 + z2 - sign(spin) * sqrt((3 - z1) * (3 + z1 + 2 * z2))
    return rms
end

"""
Computes the black hole accretion efficiency parameter from the isco.

# Parameters
- spin : black hole spin parameter ∈ (-1, 1)
"""
function compute_efficiency(spin)
    isco = compute_isco(spin)
    return 1 - sqrt(1 - 2 / (3 * isco))
end

"""
Computes the mass accretion rate to the black hole,
assuming no wind losses.

# Parameters
- M : black hole mass [in grams]
- mdot : normalized eddington rate ∈ [0,1]
- spin : black hole spin parameter ∈ (-1, 1)
"""
function compute_mass_accretion_rate(M, mdot, spin; mu_electron = 1.17)
    efficiency = compute_efficiency(spin)
    return compute_bolometric_luminosity(M, mdot, mu_electron = mu_electron) /
           (efficiency * C^2)
end

function compute_mass_accretion_rate(bh; mu_electron = 1.17)
    return compute_mass_accretion_rate(bh.M, bh.mdot, bh.spin, mu_electron = mu_electron)
end
