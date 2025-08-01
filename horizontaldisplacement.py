import numpy as np
import matplotlib.pyplot as plt
from basefunctions import *
n1=1.05
n2 = 1.00
N_photons = 10000
pore_fraction = 0.985
pore_radius = 1e-4 #1e-4 to 1e-7
foam_height = 1700e-4 #cm
foam_diameter = 900e-4 #cm
foam_radius = foam_diameter/2
foam_volume = np.pi*((foam_diameter)/2)**2*foam_height
pore_volume = 4/3*(np.pi)*pore_radius**3
N_pores = int(pore_fraction*(foam_volume)/pore_volume)

def simulate_horizontal_displacement(N_photons, foam_radius, pore_radius, pore_fraction):
    avg_pores = avg_interactions(pore_fraction, foam_height, pore_radius)
    displacements = np.zeros(N_photons)
    for i in range(N_photons):
        #N_interactions = np.random.poison(avg_interactions)
        dx_tot = 0.0
        for j in range(avg_pores):
            d = np.random.uniform(-1*pore_radius, pore_radius)
            theta_i = angle_of_incidence(d, pore_radius)
            theta_t = snells_law(theta_i, n1, n2)
            dx = lateral_displacement(theta_i, theta_t, pore_radius)
            dx_tot += dx
            theta_i = 2*theta_t - theta_i
        displacements[i] = dx_tot
    print(f'Pure attenuation would result in 0 horizontal displacements')
    print(f'The average horizontal displacmeents considering scattering is {np.mean(displacements)}')
    print(f'The percent difference of horizontal displacements caused by scattering is {np.mean(displacements)*100}')
    return displacements


displacements = simulate_horizontal_displacement(N_photons, foam_radius, pore_radius, pore_fraction)

#Plot net horizontal displacements
plt.hist(displacements *10**4, bins=100, color='skyblue', edgecolor='k')
plt.title("Simulated Photon Lateral Displacement Distribution")
plt.xlabel("Net Horizontal Displacement (μm)")
plt.ylabel("Number of Photons")
plt.tight_layout()
plt.show()


