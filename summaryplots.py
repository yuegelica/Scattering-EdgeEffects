import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as const
from mendeleev import element
from basefunctions import * 
from matplotlib.ticker import ScalarFormatter 
import matplotlib.ticker as ticker
# Constants
h = const.h
c = const.c
hc = h * c
r_e = const.physical_constants['classical electron radius'][0] * 100  # in cm
# Initial Vars
n1 = 1.05  # Refractive index inside foam
n2 = 1.00  # Refractive index of pores
foam_height = 1700e-4  # (1700 μm)
beam_width = 900e-4  # (900 μm)
rho = 2.33
energy_keV = 2.0
formula = "Si5O10Ti"
num_rays = 60  # Can be adjusted based on run time

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

def simulate_foam_metrics(pore_fraction, pore_radius, n_iterations=10):
    # Calculate attenuation
    energy_eV = energy_keV * 1000
    energy_grid, mu_values = weighted_ma(formula, rho)
    mu_at_energy = np.interp(energy_eV, energy_grid, mu_values)
    pure_trans = np.exp(-mu_at_energy * foam_height)

    # Run multiple iterations and average results
    trans_ratios = []
    displacements_list = []
    
    for iteration in range(n_iterations):
        # Initialize beam
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

        # Calculate metrics for this iteration
        path_lengths = np.sqrt(foam_height**2 + displacements**2)
        scattered_trans = np.mean(np.exp(-mu_at_energy * path_lengths))
        trans_ratios.append(scattered_trans / pure_trans)
        displacements_list.append(np.mean(np.abs(displacements)) * 1e4)  # Convert to μm
    
    # Return averaged results
    return np.mean(trans_ratios), np.mean(displacements_list)

def run_parameter_sweeps():
    # Define parameter ranges
    porosities = np.array([50, 60, 70, 80, 90, 95]) / 100  # Convert to fraction
    pore_radii = np.array([1e-7, 1e-6, 1e-5, 3e-5, 1e-4])  # cm

    # Fixed parameters
    fixed_porosity = 0.985  # 98.5%
    fixed_pore_radius = 1e-4  # 1 μm

    print("Running pore radius sweep (fixed porosity = 98.5%)...")
    radius_trans_ratios = []
    radius_displacements = []

    for radius in pore_radii:
        print(f"  Pore radius: {radius*1e4:.1f} μm")
        trans_ratio, displacement = simulate_foam_metrics(fixed_porosity, radius)
        radius_trans_ratios.append(trans_ratio)
        radius_displacements.append(displacement)

    porosity_trans_ratios = []
    porosity_displacements = []

    for porosity in porosities:
        print(f"  Porosity: {porosity*100:.1f}%")
        trans_ratio, displacement = simulate_foam_metrics(porosity, fixed_pore_radius)
        porosity_trans_ratios.append(trans_ratio)
        porosity_displacements.append(displacement)


    # Plot 1: Transmission vs Pore Radius
    plt.figure(figsize=(7, 5))
    plt.loglog(pore_radii * 1e4, radius_trans_ratios, 'bo-', linewidth=2, markersize=8)
    plt.xlabel('Pore Radius (μm)')
    plt.ylabel('Transmission Ratio (Scattered / Pure)')
    plt.title('Transmission Ratio vs Pore Radius\n(Porosity = 98.5%)')
    plt.grid(True, alpha=0.3)
    plt.xlim(left=0.001)
    plt.tight_layout()
    plt.show()


    # Plot 2: Transmission vs Porosity
    plt.figure(figsize=(7, 5))
    plt.plot(porosities * 100, porosity_trans_ratios, 'ro-', linewidth=2, markersize=8)
    plt.xlabel('Porosity (%)')
    plt.ylabel('Transmission Ratio (Scattered / Pure)')
    plt.title('Transmission Ratio vs Porosity\n(Pore Radius = 1 μm)')
    plt.grid(True, alpha=0.3)
    plt.xlim(45, 99.5)

    # Optional: Force plain number format on y-axis
    ax = plt.gca()
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f'{y:.6f}'))

    plt.tight_layout()
    plt.show()

# Plot 3: Displacement vs Pore Radius (semilogx)
    plt.figure(figsize=(7, 5))
    plt.semilogx(pore_radii * 1e4, radius_displacements, 'go-', linewidth=2, markersize=8)
    plt.xlabel('Pore Radius (μm)')
    plt.ylabel('Mean |Horizontal Displacement| (μm)')
    plt.title('Horizontal Displacement vs Pore Radius\n(Porosity = 98.5%)')
    plt.grid(True, which='both', alpha=0.3)
    plt.xlim(left=0.001)
    plt.tight_layout()
    plt.show()


    # Plot 4: Displacement vs Porosity
    plt.figure(figsize=(7, 5))
    plt.plot(porosities * 100, porosity_displacements, 'mo-', linewidth=2, markersize=8)
    plt.xlabel('Porosity (%)')
    plt.ylabel('Mean |Horizontal Displacement| (μm)')
    plt.title('Horizontal Displacement vs Porosity\n(Pore Radius = 1 μm)')
    plt.grid(True, alpha=0.3)
    plt.xlim(45, 99.5)
    plt.tight_layout()
    plt.show()

    return {
        'pore_radii': pore_radii,
        'porosities': porosities,
        'radius_trans_ratios': radius_trans_ratios,
        'radius_displacements': radius_displacements,
        'porosity_trans_ratios': porosity_trans_ratios,
        'porosity_displacements': porosity_displacements
    }


# Run the analysis
if __name__ == "__main__":
    results = run_parameter_sweeps()
