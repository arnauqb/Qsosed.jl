using Pkg
Pkg.activate("/home/arnau/code/qsosed")
using Qsosed
using PyPlot
LogNorm = matplotlib.colors.LogNorm


model = QsosedModel("./configs/config_example.yaml");

energy_range = 10 .^ range(-3, 2, length=50);

corona_flux = compute_corona_photon_flux(model, energy_range);
disk_flux = compute_disk_photon_flux(model.bh, r_min = model.warm.radius, r_max= gravity_radius(model.bh), energy_range=energy_range);

fig, ax = plt.subplots()
ax.loglog(energy_range, disk_flux, label = "disk")
ax.loglog(energy_range, corona_flux, label = "corona")
ax.set_ylim(maximum(disk_flux) / 1e2, 2 * maximum(disk_flux))
ax.legend()



radius = [10, 50, 100, 500]
fig, ax = plt.subplots()
for (i, r) in enumerate(radius)
    ph_r = compute_disk_photon_flux_at_radius(model.bh, r, energy_range)
    ax.loglog(energy_range[2:end], ph_r, label = radius[i])
end
ax.set_ylim(1e35, 1e41)
ax.legend()
