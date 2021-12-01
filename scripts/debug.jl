using Pkg
Pkg.activate("/home/arnau/code/qsosed.jl")
using Qsosed
using PyPlot
LogNorm = matplotlib.colors.LogNorm


model = QsosedModel("./configs/config_example.yaml");

energy_range = 10 .^ range(-5, 3, length = 500);
cl = compute_corona_spectral_luminosity(model.corona, model.warm, energy_range);
dl = Qsosed.compute_disk_spectral_luminosity(
    model.bh,
    energy_range,
    r_min = model.warm.radius,
    r_max = gravity_radius(model.bh),
);
wl = compute_warm_spectral_luminosity(model.corona, model.warm, energy_range);
fig, ax = plt.subplots()
dl = dl .* energy_range;
cl = cl .* energy_range;
wl = wl .* energy_range;
ax.loglog(energy_range, dl, label = "disk")
ax.loglog(energy_range, cl, label = "corona")
ax.loglog(energy_range, wl, label = "warm")
ax.set_ylim(maximum(dl) / 1e2, maximum(dl) * 1.25)
ax.set_xlim(5e-5, 5e2)
ax.legend()

total_xray_fraction(model)
