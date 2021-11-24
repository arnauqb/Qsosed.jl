using PyCall
export donthcomp

"""
This function was adapted by ADT from the subroutine donthcomp in
donthcomp.f, distributed with XSpec.
Nthcomp documentation:
https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNthcomp.html
Refs:
Zdziarski, Johnson & Magdziarz 1996, MNRAS, 283, 193,
as extended by Zycki, Done & Smith 1999, MNRAS 309, 561
Note that the subroutine has been modified so that parameter 4
is ignored, and the seed spectrum is always a blackbody.
ear: Energy vector, listing "Energy At Right" of bins (keV)
param: list of parameters; see the 5 parameters listed below.
The original fortran documentation for this subroutine is included below:
Driver for the Comptonization code solving Kompaneets equation
seed photons  -  (disk) blackbody
reflection + Fe line with smearing

Model parameters:
1: photon spectral index
2: plasma temperature in keV
3: (disk)blackbody temperature in keV
4: type of seed spectrum (0 - blackbody, 1 - diskbb)
5: redshift
"""
function donthcomp(ear, param)
    ne = length(ear)

    zfactor = 1 + param[5]

    xth, nth, spt = _thcompton(param[3] / 511, param[2] / 511, param[1])

    xninv = 511 / zfactor
    ih = 2
    xx = 1 / xninv
    while (ih < nth) && (xx > xth[ih])
        ih = ih + 1
    end
    il = ih - 1
    spp = spt[il] + (spt[ih] - spt[il]) * (xx - xth[il]) / (xth[ih] - xth[il])
    normfac = 1 / spp

    photar = zeros(ne)
    prim = zeros(ne)
    j = 1
    for i = 1:ne
        while (j <= nth) && (511.0 * xth[j] < ear[i] * zfactor)
            j = j + 1
        end
        if (j <= nth)
            if (j > 1)
                jl = j - 1
                prim[i] =
                    spt[jl] + (
                        (ear[i] / 511.0 * zfactor - xth[jl]) * (spt[jl + 1] - spt[jl]) /
                        (xth[jl + 1] - xth[jl])
                    )
            else
                prim[i] = spt[1]
            end
        end
    end
    for i = 2:ne
        photar[i] = (
            0.5 *
            (prim[i] / ear[i]^2 + prim[i - 1] / ear[i - 1]^2) *
            (ear[i] - ear[i - 1]) *
            normfac
        )
    end
    return photar
end

function _thcompton(tempbb, theta, gamma)
    """
    This function was adapted by ADT from the subroutine thcompton in
    donthcomp.f, distributed with XSpec.
    Nthcomp documentation:
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNthcomp.html
    Refs:
    Zdziarski, Johnson & Magdziarz 1996, MNRAS, 283, 193,
    as extended by Zycki, Done & Smith 1999, MNRAS 309, 561
    The original fortran documentation for this subroutine is included below:
    Thermal Comptonization; solves Kompaneets eq. with some
    relativistic corrections. See Lightman Zdziarski (1987), ApJ
    The seed spectrum is a blackbody.
    version: January 96
    #c  input parameters:
    #real * 8 tempbb,theta,gamma
    """
    #c use internally Thomson optical depth
    tautom = sqrt(2.250 + 3.0 / (theta * ((gamma + 0.50)^2 - 2.250))) - 1.50

    # Initialise arrays
    dphdot = zeros(900)
    rel = zeros(900)
    c2 = zeros(900)
    sptot = zeros(900)
    bet = zeros(900)
    x = zeros(900)

    #c JMAX  -  # OF PHOTON ENERGIES
    #c delta is the 10 - log interval of the photon array.
    delta = 0.02
    deltal = delta * log(10.0)
    xmin = 1e-4 * tempbb
    xmax = 40.0 * theta
    jmax = min(899, Int(round(log10(xmax / xmin) / delta)) + 1)

    #c X  -  ARRAY FOR PHOTON ENERGIES
    # Energy array is normalized by 511 keV, the rest energy of an electron
    for j=1:jmax+1
        x[j] = xmin * 10 ^ ((j-1) * delta)
    end

    #c compute c2(x), and rel(x) arrays
    #c c2(x) is the relativistic correction to Kompaneets equation
    #c rel(x) is the Klein - Nishina cross section divided by the
    #c Thomson crossection
    for j = 1:jmax
        w = x[j]
        #c c2 is the Cooper's coefficient calculated at w1
        #c w1 is x(j + 1 / 2) (x(i) defined up to jmax + 1)
        w1 = sqrt(x[j] * x[j + 1])
        c2[j] = (w1^4 / (1.0 + 4.60 * w1 + 1.1 * w1 * w1))
        if (w <= 0.05)
            #c use asymptotic limit for rel(x) for x less than 0.05
            rel[j] = (1.0 - 2.0 * w + 26.0 * w * w * 0.2)
        else
            z1 = (1.0 + w) / w^3
            z2 = 1.0 + 2.0 * w
            z3 = log(z2)
            z4 = 2.0 * w * (1.0 + w) / z2
            z5 = z3 / 2.0 / w
            z6 = (1.0 + 3.0 * w) / z2 / z2
            rel[j] = (0.75 * (z1 * (z4 - z3) + z5 - z6))
        end
    end

    #c the thermal emission spectrum
    jmaxth = min(900, Int(round(log10(50 * tempbb / xmin) / delta)))
    if jmaxth > jmax
        jmaxth = jmax
    end
    planck = 15.0 / (π * tempbb)^4
    for j=1:jmaxth
        dphdot[j] = planck * x[j]^2 / (exp(x[j] / tempbb)-1)
    end

    #c compute beta array, the probability of escape per Thomson time.
    #c bet evaluated for spherical geometry and nearly uniform sources.
    #c Between x = 0.1 and 1.0, a function flz modifies beta to allow
    #c the increasingly large energy change per scattering to gradually
    #c eliminate spatial diffusion
    jnr = Int(round(log10(0.10 / xmin) / delta + 1))
    jnr = min(jnr, jmax - 1)
    jrel = Int(round(log10(1 / xmin) / delta + 1))
    jrel = min(jrel, jmax)
    xnr = x[jnr]
    xr = x[jrel]
    for j = 1:(jnr - 1)
        taukn = tautom * rel[j]
        bet[j] = 1.0 / tautom / (1.0 + taukn / 3.0)
    end
    for j = jnr:jrel
        taukn = tautom * rel[j]
        arg = (x[j] - xnr) / (xr - xnr)
        flz = 1 - arg
        bet[j] = 1.0 / tautom / (1.0 + taukn / 3.0 * flz)
    end
    for j = jrel+1:jmax
        bet[j] = 1.0 / tautom
    end

    dphesc = _thermlc(tautom, theta, deltal, x, jmax, dphdot, bet, c2)

    #c     the spectrum in E F_E
    for j = 1:(jmax - 1)
        sptot[j] = dphesc[j] * x[j]^2
    end

    return x, jmax, sptot
end

function _thermlc(tautom, theta, deltal, x, jmax, dphdot, bet, c2)
    """
    This function was adapted by ADT from the subroutine thermlc in
    donthcomp.f, distributed with XSpec.
    Nthcomp documentation:
    https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSmodelNthcomp.html
    Refs:
    Zdziarski, Johnson & Magdziarz 1996, MNRAS, 283, 193,
    as extended by Zycki, Done & Smith 1999, MNRAS 309, 561
    The original fortran documentation for this subroutine is included below:
    This program computes the effects of Comptonization by
    nonrelativistic thermal electrons in a sphere including escape, and
    relativistic corrections up to photon energies of 1 MeV.
    the dimensionless photon energy is x = hv / (m * c * c)
    The input parameters and functions are:
    dphdot(x), the photon production rate
    tautom, the Thomson scattering depth
    theta, the temperature in units of m*c*c
    c2(x), and bet(x), the coefficients in the K - equation and the
      probability of photon escape per Thomson time, respectively,
      including Klein - Nishina corrections
    The output parameters and functions are:
    dphesc(x), the escaping photon density
    """
    dphesc = zeros(900)  # Initialise the output
    a = zeros(900)
    b = zeros(900)
    c = zeros(900)
    d = zeros(900)
    alp = zeros(900)
    u = zeros(900)
    g = zeros(900)
    gam = zeros(900)

    #c u(x) is the dimensionless photon occupation number
    c20 = tautom / deltal

    #c determine u
    #c define coefficients going into equation
    #c a(j) * u(j + 1) + b(j) * u(j) + c(j) * u(j - 1) = d(j)
    for j = 2:(jmax - 1)
        w1 = sqrt(x[j] * x[j + 1])
        w2 = sqrt(x[j - 1] * x[j])
        #c  w1 is x(j + 1 / 2)
        #c  w2 is x(j - 1 / 2)
        a[j] = -c20 * c2[j] * (theta / deltal / w1 + 0.5)
        t1 = -c20 * c2[j] * (0.5 - theta / deltal / w1)
        t2 = c20 * c2[j - 1] * (theta / deltal / w2 + 0.5)
        t3 = x[j]^3 * (tautom * bet[j])
        b[j] = t1 + t2 + t3
        c[j] = c20 * c2[j - 1] * (0.5 - theta / deltal / w2)
        d[j] = x[j] * dphdot[j]
    end

    #c define constants going into boundary terms
    #c u(1) = aa * u(2) (zero flux at lowest energy)
    #c u(jx2) given from region 2 above
    x32 = sqrt(x[1] * x[2])
    aa = (theta / deltal / x32 + 0.5) / (theta / deltal / x32 - 0.5)

    #c zero flux at the highest energy
    u[jmax] = 0.0

    #c invert tridiagonal matrix
    alp[2] = b[2] + c[2] * aa
    gam[2] = a[2] / alp[2]
    for j = 3:(jmax - 1)
        alp[j] = b[j] - c[j] * gam[j - 1]
        gam[j] = a[j] / alp[j]
    end
    g[2] = d[2] / alp[2]
    for j = 3:(jmax - 2)
        g[j] = (d[j] - c[j] * g[j - 1]) / alp[j]
    end
    g[jmax - 1] =
        (d[jmax - 1] - a[jmax - 1] * u[jmax] - c[jmax-1] * g[jmax - 2]) /
        alp[jmax - 1]
    u[jmax - 1] = g[jmax - 1]
    for j = 3:jmax-1
        jj = jmax + 1 - j
        u[jj] = g[jj] - gam[jj] * u[jj + 1]
    end
    u[1] = aa * u[2]
    #c compute new value of dph(x) and new value of dphesc(x)
    dphesc[1:jmax] = @. x[1:jmax] * x[1:jmax] * u[1:jmax] * bet[1:jmax] * tautom

    return dphesc
end
