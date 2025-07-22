import numpy as np

from tools.parse_sim import parse_sim
from swiftgalaxy import SWIFTGalaxy, SOAP

args, config = parse_sim()

sg = SWIFTGalaxy(
    config['ColibrePaths']['virtualSnapshotPath'],
    SOAP(config['ColibrePaths']['haloCataloguePath'],
        soap_index=184151)
)

gas_coords = sg.gas.coordinates.to_physical().to('kpc')  # Convert coordinates to physical units in kpc
print(gas_coords.units)
dust_r = np.sqrt(np.sum(gas_coords**2, axis=1))  # Calculate the radius in pc
dust_m = (sg.gas.total_dust_mass_fractions * sg.gas.masses).to_physical().to('Msun')  # Dust mass fraction in physical units
dustMasses_sorted = dust_m[np.argsort(dust_r)]

idx_halfmass = np.min(np.argwhere((np.cumsum(dustMasses_sorted) / np.sum(dustMasses_sorted)) >= 0.5))

dustHalfMassRadius = np.sort(dust_r)[idx_halfmass]
dustHalfMass = (np.sum(dust_m) / 2.)


SigmaDust = dustHalfMass / (np.pi * dustHalfMassRadius**2) # In solar masses / kpc^2

maxDustFraction = np.clip(10**(-0.5 - np.log10(SigmaDust)), a_min = 10**(-6.5), a_max = 10**(-4.5))

print(f"Max dust fraction: {maxDustFraction:.6e} for galaxy with SOAP index {sg.halo_catalogue.soap_index}")

# def max_dust_fraction(sg):
#     """Compute the maximum dust fraction for each galaxy in the SWIFTGalaxy object."""
#     gas_coords = sg.gas.coordinates.to_physical().to('pc')
#     gas_r = np.sqrt(np.sum(gas_coords**2, axis=1))  # Calculate the radius in pc
    
#     if gas_r.size == 0:
#         return 10**(-8.5)
    
#     dust_masses = (sg.gas.total_dust_mass_fractions * sg.gas.masses).to_physical().to('Msun')  # Dust mass fraction in physical units
#     total_dust_halfmass = np.sum(dust_masses) / 2.0   # Total dust mass in Msun
#     # Sort the dust masses by radius
#     sorted_indices = np.argsort(gas_r)
#     sorted_dust_masses = dust_masses[sorted_indices]
#     sorted_gas_r = gas_r[sorted_indices]

#     # Find the index of the half-mass radius
#     cumulative_dust_mass = np.cumsum(sorted_dust_masses)
#     half_mass_index = np.searchsorted(cumulative_dust_mass, total_dust_halfmass)

#     dust_half_mass_radius = sorted_gas_r[half_mass_index]  # In kpc

#     print(total_dust_halfmass, dust_half_mass_radius)
#     sigma_dust = total_dust_halfmass / (np.pi * dust_half_mass_radius**2) # In solar masses / pc^2
    
#     max_dust_fraction = np.clip(10**(-0.5 - np.log10(sigma_dust)), a_min=10**(-6.5), a_max=10**(-4.5))
    
#     print(f"Calculating maximum dust fraction for galaxy {sg.halo_catalogue.soap_index} with max dust fraction: {max_dust_fraction:.6e}")
#     return max_dust_fraction

print("With correct method: ", maxDustFraction)
# print("With function: ", max_dust_fraction(sg))
