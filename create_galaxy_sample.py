"""For each snap we create a sample of galaxies which will be used for the SKIRT run.
These galaxies will go through certain selection criteria and will output a sample file

This sample file will contain the halo IDs, stellar masses, and max dust fractions of each galaxy"""
import os
from swiftsimio import load as load_snapshot
from swiftgalaxy.iterator import SWIFTGalaxies
from swiftgalaxy.halo_catalogues import SOAP

from unyt import unyt_array
import numpy as np 
from tools.parse_sim import parse_sim

args, config = parse_sim()

# Now we set the selection criteria for the mass. Make edits in the SKIRTconfig.yaml file.
mass_lb, mass_ub = unyt_array(list(config['SelectionCriteria']['massRange'].values()), units='Msun', dtype=float) # add dtype float to avoid unyt_array error

def max_dust_fraction(sg):
    """Compute the maximum dust fraction for each galaxy in the SWIFTGalaxy object."""
    gas_coords = sg.gas.coordinates.to_physical().to('kpc')
    gas_r = np.sqrt(np.sum(gas_coords**2, axis=1))  # Calculate the radius in pc
    
    if gas_r.size == 0:
        return 0. # If there are no gas particles, return 0 dust fraction
    elif gas_r.size == 1:
        return 10**(-4.5) # If there is only one gas particle, return the maximum dust fraction
    else:
        dust_masses = (sg.gas.total_dust_mass_fractions * sg.gas.masses).to_physical().to('Msun')  # Dust mass fraction in physical units
        total_dust_halfmass = np.sum(dust_masses) / 2.0   # Total dust mass in Msun
        
        # Sort the dust masses by radius
        sorted_indices = np.argsort(gas_r)
        sorted_dust_masses = dust_masses[sorted_indices]
        sorted_gas_r = gas_r[sorted_indices]

        # Find the index of the half-mass radius
        cumulative_dust_mass = np.cumsum(sorted_dust_masses)
        half_mass_index = np.searchsorted(cumulative_dust_mass, total_dust_halfmass)
        dust_half_mass_radius = sorted_gas_r[half_mass_index]  # In kpc

        sigma_dust = total_dust_halfmass / (np.pi * dust_half_mass_radius**2) # In solar masses / pc^2
        max_dust_fraction = np.clip(10**(-0.5 - np.log10(sigma_dust)), a_min=10**(-6.5), a_max=10**(-4.5))
        
        print(f"Calculating maximum dust fraction for galaxy {sg.halo_catalogue.soap_index} with max dust fraction: {max_dust_fraction:.6e}")
        return max_dust_fraction

# Load the simulation paths
virtual_snap_file = config['ColibrePaths']['virtualSnapshotPath']
catalogue_file = config['ColibrePaths']['haloCataloguePath']
catalogue = load_snapshot(catalogue_file)

# Establish the path where the SKIRT simulation is stored
skirt_dir = config['SKIRTPaths']['skirtRunsPath']
if not os.path.exists(skirt_dir):
    os.makedirs(skirt_dir)
    print(f"Created directory: {skirt_dir}")
else:
    print(f"Directory already exists: {skirt_dir}")

# Create the information for the sample file
halo_IDs = catalogue.input_halos.halo_catalogue_index.value
halo_track_IDs = catalogue.input_halos_hbtplus.track_id.value
halo_stellar_mass = unyt_array(catalogue.exclusive_sphere_50kpc.stellar_mass.to_physical())

mass_selection = (halo_stellar_mass > mass_lb) * (halo_stellar_mass < mass_ub)
halo_indices = np.where(mass_selection)[0]  # Get the indices of the selected halos

# Calculate the maximum dust fraction for each galaxy
sgs = SWIFTGalaxies(virtual_snap_file, halo_catalogue=SOAP(catalogue_file, soap_index=halo_indices), preload={'gas.coordinates', 'gas.smoothing_lengths', 'gas.masses', 'gas.metal_mass_fractions', 'gas.star_formation_rates', 'gas.averaged_star_formation_rates','gas.total_dust_mass_fractions','gas.dust_mass_fractions.GraphiteLarge', 'gas.dust_mass_fractions.MgSilicatesLarge', 'gas.dust_mass_fractions.FeSilicatesLarge','gas.dust_mass_fractions.GraphiteSmall', 'gas.dust_mass_fractions.MgSilicatesSmall', 'gas.dust_mass_fractions.FeSilicatesSmall'})
print(f"Loaded {mass_selection.sum()} galaxies with stellar mass between {mass_lb:.6e} and {mass_ub:.6e} Msun.")
max_dust_fractions = sgs.map(max_dust_fraction)

# Create the sample file with the selected galaxies

sample_info = f"""Simulation: {args.sim_name}
Snapshot number: {args.snap_nr}
Mass selection: {mass_lb:.6e} < Stellar Mass < {mass_ub:.6e}
Number of galaxies: {mass_selection.sum()}"""

sample_header = """Column 1: Halo IDs
Column 2: Halo Track IDs
Column 3: Stellar Mass (Msun)
Column 4: Maximum Dust Fraction (1)"""

sample_file = np.column_stack((halo_IDs, halo_track_IDs, halo_stellar_mass.to('Msun').value))[mass_selection, :]
sample_file = np.column_stack((sample_file, max_dust_fractions))
sample_file_path = os.path.join(skirt_dir, f"sample_{args.snap_nr}.txt")
np.savetxt(sample_file_path, sample_file, header=sample_header, fmt='%d %d %.6e %.6e', delimiter=' ', comments='#')
np.savetxt(os.path.join(skirt_dir, f"halo_indices_{args.snap_nr}.txt"), halo_indices, fmt='%d', header='Halo indices of selected galaxies')

# Create a sample info file for the user and write the information to it
with open(os.path.join(skirt_dir, f"sample_info_{args.snap_nr}.txt"), 'w') as f:
    f.write(sample_info)