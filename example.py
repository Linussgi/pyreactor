import numpy as np
from matplotlib import pyplot as plt

from pyreactor.reaction import Reaction
from pyreactor.component import Component

# Example reaction using the Haber process

# Define shomate coefficients
N2_COEFFS = [19.50583, 19.88705, -8.598535, 1.369784, 0.527601]
H2_COEFFS = [33.066178, -11.363417, 11.432816, -2.772874, -0.158558]
NH3_COEFFS = [19.99563, 49.77119, -15.37599, 1.921168, 0.189174]

# Define all reaction components
ammonia = Component("NH3", 2, NH3_COEFFS)
hydrogen = Component("H2", 3, H2_COEFFS)
nitrogen = Component("N2", 1, N2_COEFFS)

# Create reaction object
haber_process = Reaction([nitrogen, hydrogen], [ammonia], -92000, -199, "ideal", 500)

# Define temperature and pressure to model reaction over
pressures = np.linspace(10, 300, 100)
temps = np.linspace(300, 2000, 100)

conv_solutions = []
for temp in temps:
    haber_process.reac_temp = temp

    k_eq = haber_process.calculate_rxn_k()

    pressure_sol = []
    for pressure in pressures:
        chi = haber_process.calculate_conversion(pressure, [100, 300], [0], k_eq)

        n2_conv = (nitrogen.order * chi) / 100

        pressure_sol.append(n2_conv)

    conv_solutions.append(pressure_sol)

# Plot results
fig, ax = plt.subplots(figsize=[8, 6])

pressure_mesh, temp_mesh = np.meshgrid(pressures, temps)
conversion_array = np.array(conv_solutions)

heatmap = ax.pcolormesh(temp_mesh, pressure_mesh, conversion_array, cmap="magma", shading="auto")

cbar = plt.colorbar(heatmap)
cbar.set_label("N2 Conversion")
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Pressure (atm)")

# plt.savefig("conversion_heatmap.svg", format="svg", bbox_inches="tight")
plt.show()