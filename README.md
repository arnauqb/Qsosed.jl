[![Build Status](https://github.com/arnauqb/qsosed.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/arnauqb/qsosed.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/arnauqb/qsosed.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/arnauqb/qsosed.jl)
[![Documentation](https://github.com/arnauqb/qsosed.jl/actions/workflows/docs.yaml/badge.svg)](https://arnauqb.github.io/Qsosed.jl)

# Qsosed.jl

This Julia package handles various calculations involving the accretion physics of AGNs. In particular, it implements the [qsosed]("https://github.com/HEASARC/xspec_localmodels/tree/master/agnsed") model of Xspec, explained in [Kubota & Done (2018)]("https://arxiv.org/abs/1804.00171") to create the flux energy distribution in the UV/X-Ray band of an AGN.

The SED model has three characteristic regions: the outer standard
disc region; the warm Comptonising region; and the inner hot
Comptonising region.

For the warm Comptonising region, this model adopts the passive disc
scenario tested by Petrucci et al. 2018
(https://ui.adsabs.harvard.edu//#abs/2018A&A...611A..59P/abstract). Here,
the flow is assumed to be completely radially stratified, emitting as
a standard disc blackbody from Rout to Rwarm, as warm Comptonisation
from Rwarm to Rhot and then makes a transition to the hard X-ray
emitting hot Comptonisation component from Rhot to RISCO. The warm
Comptonisation component is optically thick, so is associated with
material in the disc. Nonetheless, the energy does not thermalise to
even a modified blackbody, perhaps indicating that significant
dissipation takes place within the vertical structure of the disc,
rather than being predominantly released in the midplane.

At a radius below Rhot, the energy is emitted in the hot
Comptonisation component. This has much lower optical depth, so it is
not the disc itself. In the model, the albedo is fixed at a = 0.3, and
the seed photon temperature for the hot Comptonisation component is
calculated internally. In contrast to optxagnf, this model does not
take the color temperature correction into account.
