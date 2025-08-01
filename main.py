import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as const
from mendeleev import element
from basefunctions import * 

# Constants
h = const.h
c = const.c
hc = h * c
r_e = const.physical_constants['classical electron radius'][0] * 100  # in cm

# Initial Variables
n1 = 1.05  # Refractive index inside foam
n2 = 1.00  # Refractive index of pores
pore_fraction = 0.985
pore_radius = 1e-4  # cm (1 μm)
foam_height = 1700e-4  # cm (1700 μm)
beam_width = 900e-4  # cm (900 μm)
ccd_pixel_size = 13.5e-4  # cm (13.5 μm)
rho = 2.33  # Density (g/cm³)


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


def run_simulation(energy_keV=2.0, formula="Si5O10Ti", pore_radius=pore_radius, pore_fraction=pore_fraction,num_rays=1000):
    # Convert energy and calculate attenuation
    energy_eV = energy_keV * 1000
    energy_grid, mu_values = weighted_ma(formula, rho)
    mu_at_energy = np.interp(energy_eV, energy_grid, mu_values)
    pure_trans = np.exp(-mu_at_energy * foam_height)

    # Initialize rays
    x_pos = np.linspace(-beam_width/2, beam_width/2, num_rays, endpoint=False)
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
    
    # Calculate positions
    beam_half_width = beam_width / 2 * 1e4  # μm (450)
    final_positions = x_pos * 1e4 + displacements * 1e4
    
    return {
        'x_pos': x_pos,'displacements': displacements,'transmissions': transmissions,'pure_trans': pure_trans,'final_positions': final_positions,'beam_half_width': beam_half_width,'avg_displacement': np.mean(np.abs(displacements)) * 1e4,'trans_reduction': 100 * (1 - np.mean(transmissions) / pure_trans),'pore_radius': pore_radius,'pore_fraction': pore_fraction
    }

results = run_simulation(num_rays=10000)  

# Plot 1: Transmission Difference
plt.figure(figsize=(10, 6))
trans_diff = results['pure_trans'] - results['transmissions']
plt.hist(trans_diff, bins=50, color='purple')
plt.axvline(0, color='k', linestyle='--')
plt.title('Transmission Difference (Pure - Scattered)')
plt.xlabel('Δ Transmission')
plt.ylabel('Counts')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Plot 2: Zoomed Left Edge
plt.figure(figsize=(10, 6))
left_region_start = -results['beam_half_width'] - 30
left_region_end = left_region_start + 100
pixel_edges_zoom = np.arange(left_region_start, left_region_end + ccd_pixel_size*1e4, ccd_pixel_size*1e4)
pixel_centers_zoom = (pixel_edges_zoom[1:] + pixel_edges_zoom[:-1])/2

pure_binned_zoom, _ = np.histogram(results['x_pos']*1e4, bins=pixel_edges_zoom, weights=np.full(len(results['x_pos']), results['pure_trans']))
scat_binned_zoom, _ = np.histogram(results['final_positions'], bins=pixel_edges_zoom, weights=results['transmissions'])

plt.bar(pixel_centers_zoom, pure_binned_zoom, width=13.5, color='blue', alpha=0.5, label='Pure')
plt.bar(pixel_centers_zoom, scat_binned_zoom, width=13.5, color='orange', alpha=0.5, label='Scattered')
for edge in pixel_edges_zoom:
    plt.axvline(edge, color='gray', alpha=0.3, linestyle=':', linewidth=0.5)
plt.axvline(-results['beam_half_width'], color='red', linestyle='--', alpha=0.7, label='Beam Edge')
plt.title('Zoomed Left Edge Comparison')
plt.xlabel('Position (μm)')
plt.ylabel('Summed Intensity')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Plot 3: Zoomed Right Edge
plt.figure(figsize=(10, 6))
right_region_start = results['beam_half_width'] - 100
right_region_end = results['beam_half_width'] + 30
pixel_edges_zoom = np.arange(right_region_start, right_region_end + ccd_pixel_size*1e4, ccd_pixel_size*1e4)
pixel_centers_zoom = (pixel_edges_zoom[1:] + pixel_edges_zoom[:-1])/2

pure_binned_zoom, _ = np.histogram(results['x_pos']*1e4, bins=pixel_edges_zoom, weights=np.full(len(results['x_pos']), results['pure_trans']))
scat_binned_zoom, _ = np.histogram(results['final_positions'], bins=pixel_edges_zoom, weights=results['transmissions'])

plt.bar(pixel_centers_zoom, pure_binned_zoom, width=13.5, color='blue', alpha=0.5, label='Pure')
plt.bar(pixel_centers_zoom, scat_binned_zoom, width=13.5, color='orange', alpha=0.5, label='Scattered')
for edge in pixel_edges_zoom:
    plt.axvline(edge, color='gray', alpha=0.3, linestyle=':', linewidth=0.5)
plt.axvline(results['beam_half_width'], color='red', linestyle='--', alpha=0.7, label='Beam Edge')
plt.title('Zoomed Right Edge Comparison')
plt.xlabel('Position (μm)')
plt.ylabel('Summed Intensity')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()

# Plot 4: Overlaid CCD 
plt.figure(figsize=(12, 6))
# Create full CCD bins
pixel_edges_full = np.arange(-results['beam_half_width'] - 100, results['beam_half_width'] + 100 + ccd_pixel_size*1e4, ccd_pixel_size*1e4)
pixel_centers_full = (pixel_edges_full[1:] + pixel_edges_full[:-1])/2

pure_binned_full, _ = np.histogram(results['x_pos']*1e4, bins=pixel_edges_full, weights=np.full(len(results['x_pos']), results['pure_trans']))
scat_binned_full, _ = np.histogram(results['final_positions'], bins=pixel_edges_full, weights=results['transmissions'])

plt.bar(pixel_centers_full, pure_binned_full, width=13.5, color='blue', alpha=0.5, label='Pure')
plt.bar(pixel_centers_full, scat_binned_full, width=13.5, color='orange', alpha=0.5, label='Scattered')
for edge in pixel_edges_full:
    plt.axvline(edge, color='gray', alpha=0.3, linestyle=':', linewidth=0.5)
plt.axvline(-results['beam_half_width'], color='red', linestyle='--', alpha=0.7, label='Beam Edge')
plt.axvline(results['beam_half_width'], color='red', linestyle='--', alpha=0.7)
plt.title('CCD Images Overlay (Pure vs Scattered)')
plt.xlabel('Position (μm)')
plt.ylabel('Summed Intensity')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
