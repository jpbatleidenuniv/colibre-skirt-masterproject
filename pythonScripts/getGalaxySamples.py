"""
Script to select COLIBRE halos and store some global information.
Created by Andrea Gebek on 25.3.2025.

Edit: In our case, snapshots are located on the Sterrenwacht (STRW) server of Leiden Observatory.
A typical path to the simulation data is: /net/hypernova/data2/COLIBRE/L0100N1504/Thermal_non_equilibrium/SOAP/halo_properties_0127.hdf5.
Associated edits to the filepaths in the SKIRT_parameters.yml file are required.
"""
import sys
import unyt
from swiftsimio import load as load_snapshot
import numpy as np
import yaml
import argparse
import os

# Set simName
parser = argparse.ArgumentParser(
    description="Select COLIBRE halos and store some global information."
)

parser.add_argument(
    "BoxSize",
    type=int,
    help="Boxsize of the simulation in Mpc.",
)

parser.add_argument( # required since simulation formats are different e.g. L0100N1504, L0025N0752
    "NumParticles",
    type=int,
    help="Number of particles in each dimension of the simulation. Similar meaning to resolution.",
    )

# parser.add_argument(
#     "Resolution",
#     type=int,
#     help="Particle mass resolution of the simulation in log10(M/Msun).",
# )

parser.add_argument(
    "--snaps",
    type=int,
    required=True,
    nargs='+',
    help="<Required> Snapshot number(s).",
)

parser.add_argument(
    "--mode",
    type=str,
    default="Thermal_non_equilibrium", # Thermal AGN feedback with non-equilibrium chemistry
    help="Simulation mode (default: Thermal).",
)

args = parser.parse_args()

sim = 'L{:04.0f}N{:04.0f}'.format(args.BoxSize, args.NumParticles) # L0025N0752, L0100N1504, L0250N0752, L0250N1504, L0500N1504
simName = sim + '/' + args.mode # adds the simulation mode to the simulation name, e.g. L0100N1504/Thermal_non_equilibrium

# Define filepaths from parameter file
dir_path = os.path.dirname(os.path.realpath(__file__))
with open(f'{dir_path}/../SKIRT_parameters.yml','r') as stream:
    params = yaml.safe_load(stream)

simPath = params['ColibreFilepaths']['simPath'].format(simName=simName)
sampleFolder = params['ColibreFilepaths']['sampleFolder'].format(simPath=simPath)

# Check if the sample folder exists
if not os.path.exists(simPath):
    print(f"Error: The specified simulation path does not exist: {simPath}")
    sys.exit(1)

if not os.path.exists(sampleFolder):
    print(f"Error: The specified sample folder does not exist: {sampleFolder}")
    sys.exit(1)

header = 'Column 1: Halo ID\n' + \
         'Column 2: Stellar mass (Msun)\n' + \
         'Column 3: Stellar half-mass radius (kpc)\n'

for snap in args.snaps:
    
    catalogue_file = params['ColibreFilepaths']['catalogueFile'].format(simPath=simPath,snap_nr=snap)
    catalogue = load_snapshot(catalogue_file)
   
    halo_IDs = catalogue.input_halos.halo_catalogue_index.value

    Mstar = unyt.unyt_array(catalogue.bound_subhalo.stellar_mass.to_physical())
    Rstar = unyt.unyt_array(catalogue.bound_subhalo.half_mass_radius_stars.to_physical())

    SEL = (Mstar >= unyt.unyt_quantity(float(params['SelectionCriteria']['minStellarMass']), 'Msun')) * (Mstar <= unyt.unyt_quantity(float(params['SelectionCriteria']['maxStellarMass']), 'Msun')) # Simple stellar mass selection. Replace this with 
    # your selection criteria.

    print(len(SEL[SEL]), 'galaxies selected in snapshot', snap)

    sample_file = np.vstack((halo_IDs, Mstar.to('Msun').value, Rstar.to('kpc').value)).T[SEL, :]

    np.savetxt(sampleFolder + 'sample_' + str(snap) + '.txt', sample_file, fmt = ['%d', '%.6e', '%.4f'], header = header)