import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as const
from mendeleev import element
from basefunctions import *  

# ============== PHYSICAL CONSTANTS ==============
h = const.h
c = const.c
hc = h * c
r_e = const.physical_constants['classical electron radius'][0] * 100  # in cm

# ============== SIMULATION PARAMETERS ==============
n1 = 1.05  # Refractive index inside foam
n2 = 1.00  # Refractive index of pores
pore_fraction = 0.985
pore_radius = 1e-4  # cm (1 μm)
foam_height = 1700e-4  # cm (1700 μm)
foam_diameter = 900e-4  # cm (900 μm)
beam_width = 900e-4  # cm (900 μm)
ccd_pixel_size = 13.5e-4  # cm (13.5 μm)
rho = 2.33

# ============== PHYSICS FUNCTIONS ==============

def weighted_ma(formula, rho):
    elements, counts = parse_compounds(formula)
    mu_rho_list = []
    energy_common = None

    for elem in elements:
        filename = f"{elem.lower()}.txt"
        energy, f1, f2, lambda1 = load_file(filename)
        mu = attenuation(f2, f1, lambda1, elem)
        mu_rho = mass_attenuation(mu, 1.0)

        if energy_common is None:
            energy_common = energy
            mu_rho_list.append(mu_rho)
        else:
            mu_rho_interp = np.interp(energy_common, energy, mu_rho)
            mu_rho_list.append(mu_rho_interp)

    mu_rho_array = np.array(mu_rho_list)
    molar_masses = np.array([element(e).mass for e in elements])
    counts = np.array(counts)
    weight_fractions = (molar_masses * counts) / np.sum(molar_masses * counts)
    weighted_mu_rho = np.sum(weight_fractions[:, None] * mu_rho_array, axis=0)
    return energy_common, weighted_mu_rho * rho


def simulate_foam(energy_keV=2.0, formula="Si5O10Ti"):
    # Convert energy and calculate attenuation
    energy_eV = energy_keV * 1000
    energy_grid, mu_values = weighted_ma(formula, rho)
    mu_at_energy = np.interp(energy_eV, energy_grid, mu_values)
    pure_trans = np.exp(-mu_at_energy * foam_height)

    # Initialize rays
    ray_spacing = 0.1 * pore_radius
    x_pos = np.linspace(-beam_width/2, beam_width/2, int(beam_width/ray_spacing), endpoint=False)
    num_rays = len(x_pos)
    print(f'The number of rays is {num_rays}')
    avg_pores = avg_interactions(pore_fraction, foam_height, pore_radius)

    # Simulate scattering
    displacements = np.zeros(num_rays)
    theta_i = np.zeros(num_rays)
    
    for _ in range(avg_pores):
        d = np.random.uniform(-pore_radius, pore_radius, size=num_rays)
        theta_i = angle_of_incidence(d, pore_radius) + theta_i
        theta_t = snells_law(theta_i, n1, n2)
        dx = lateral_displacement(theta_i, theta_t, pore_radius)
        displacements += dx
        theta_i = 2 * theta_t - theta_i

    # Calculate transmissions
    path_lengths = np.sqrt(foam_height**2 + displacements**2)
    transmissions = np.exp(-mu_at_energy * path_lengths)
    print(transmissions)
    trans_diff = transmissions - pure_trans

    beam_half_width = beam_width / 2 * 1e4  # 450 microns
    max_displacement = np.max(np.abs(displacements)) * 1e4  # in microns
    plot_range = max(beam_half_width * 1.5, max_displacement * 1.5) 
    pixel_edges = np.arange(-plot_range, plot_range, 13.5)  #bin 13.5 microns
    pixel_centers = (pixel_edges[1:] + pixel_edges[:-1])/2

    plt.figure(figsize=(15, 12))

    # 1. Input beam profile
    plt.subplot(321)
    plt.bar(x_pos*1e4, np.ones(num_rays), width=ray_spacing*1e4, 
            color='blue', edgecolor='black', alpha=0.7)
    plt.title(f'1. Input Beam Profile ({ray_spacing*1e4:.1f}μm spacing)')
    plt.xlabel('Position (μm)')
    plt.ylabel('Intensity')
    plt.xlim(-plot_range, plot_range)
    plt.axvline(-beam_half_width, color='red', linestyle='--', alpha=0.5)
    plt.axvline(beam_half_width, color='red', linestyle='--', alpha=0.5)
    plt.text(0, 0.5, f'Beam Diameter: {beam_width*1e4:.0f} μm', 
            ha='center', va='center', backgroundcolor='white')

    # 2. Pure attenuation CCD
    plt.subplot(322)
    pure_binned = np.histogram(x_pos*1e4, bins=pixel_edges, 
                            weights=np.full(num_rays, pure_trans))[0]
    plt.bar(pixel_centers, pure_binned, width=13.5, color='blue', edgecolor='white')
    for edge in pixel_edges:
        plt.axvline(edge, color='gray', alpha=0.3, linestyle=':', linewidth=0.5)
    plt.title('2. Pure Attenuation CCD')
    plt.xlabel('Position (μm)')
    plt.ylabel('Intensity')
    plt.xlim(-plot_range, plot_range)
    plt.axvline(-beam_half_width, color='red', linestyle='--', alpha=0.5)
    plt.axvline(beam_half_width, color='red', linestyle='--', alpha=0.5)

    # 3. Displacement histogram 
    plt.subplot(323)
    plt.hist(displacements*1e4, bins=50, color='green')
    plt.title('3. Horizontal Displacements')
    plt.xlabel('Displacement (μm)')
    plt.ylabel('Counts')
    plt.grid(True, alpha=0.3)

    # 4. Transmission distribution 
    plt.subplot(324)
    plt.hist(transmissions, bins=50, color='orange')
    plt.axvline(pure_trans, color='r', linestyle='--', label=f'Pure: {pure_trans:.8f}')
    plt.axvline(np.mean(transmissions), color='b', linestyle='-', 
                label=f'Mean: {np.mean(transmissions):.8f}')
    plt.title('4. Transmission Distribution')
    plt.xlabel('Transmission')
    plt.ylabel('Counts')
    plt.gca().ticklabel_format(axis='x', style='plain', useOffset=False)  # Fixed formatting
    plt.legend()

    # 5. Scattered CCD (with reduced intensity)
    plt.subplot(325)
    final_positions = x_pos*1e4 + displacements*1e4
    scat_binned = np.histogram(final_positions, bins=pixel_edges, weights=transmissions)[0]
    plt.bar(pixel_centers, scat_binned, width=13.5, color='orange', edgecolor='white', alpha=0.8)
    for edge in pixel_edges:
        plt.axvline(edge, color='gray', alpha=0.3, linestyle=':', linewidth=0.5)
    plt.title('5. Scattered CCD Image')
    plt.xlabel('Position (μm)')
    plt.ylabel('Intensity')
    plt.xlim(-plot_range, plot_range)
    plt.axvline(-beam_half_width, color='red', linestyle='--', alpha=0.5)
    plt.axvline(beam_half_width, color='red', linestyle='--', alpha=0.5)

    # 6. Final position histogram 
    plt.subplot(326)
    plt.hist(final_positions, bins=50, color='purple')
    plt.title('6. Final Positions of Rays')
    plt.xlabel('Position (μm)')
    plt.ylabel('Counts')
    plt.grid(True, alpha=0.3)
    plt.axvline(-beam_half_width, color='red', linestyle='--', alpha=0.5)
    plt.axvline(beam_half_width, color='red', linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.show()

# Run simulation
simulate_foam(energy_keV=2.0, formula="Si5O10Ti")

