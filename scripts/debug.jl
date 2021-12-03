using Pkg
Pkg.activate("/home/arnau/code/qsosed.jl")
using Qsosed
using PyPlot, YAML, Printf
LogNorm = matplotlib.colors.LogNorm


fpath = "./configs/config_example.yaml"
config = YAML.load_file(fpath, dicttype=Dict{Symbol, Any})
model = QsosedModel(config);


M_range = collect(10 .^ (6:10));
mdot_range = 10 .^ range(log10(0.025), log10(0.9), length=20);
xray_fractions = zeros(length(M_range), length(mdot_range))
for (i, M) in enumerate(M_range)
    for (j, mdot) in enumerate(mdot_range)
        conf = copy(config)
        conf[:M] = M
        conf[:mdot] = mdot
        model = QsosedModel(conf)
        fx = total_xray_fraction(model)
        xray_fractions[i, j] = fx
    end
end

fig, ax = plt.subplots()
cm = ax.pcolormesh(log10.(M_range), log10.(mdot_range), xray_fractions')
cbar = plt.colorbar(cm, ax=ax)
cbar.set_label("X-ray fraction", rotation=-90, labelpad=15)
ax.set_xticks(6:10)
ax.set_xticklabels([L"10^{%$i}" for i in 6:10])
ax.set_yticks(log10.(mdot_range))
ax.set_yticklabels([@sprintf "%.3f" mdot for mdot in mdot_range])


total_xray_fraction(model)

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
ax.axvline(0.06, color = "black", alpha=0.5)
ax.axvline(1, color = "black", alpha=0.5)
ax.legend()

