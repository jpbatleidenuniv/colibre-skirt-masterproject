"""
Script to select COLIBRE halos and store some global information.
Created by Andrea Gebek on 25.3.2025.
"""

import unyt
from swiftsimio import load as load_snapshot
import numpy as np

dataPath = '/cosma8/data/dp004/colibre/Runs/'
sampleFolder = '/cosma/home/do019/dc-gebe1/' # Folder where the galaxy sample .txt files are stored.

simL = 25 # Box length in Mpc
simR = 6 # Mass resolution in log10(M/Msun)
simName = 'Thermal' # Thermal AGN feedback with non-equilibrium chemistry
snapList = [56, 123] # List of snapshots

sim = 'L{:03.0f}_m{:01.0f}'.format(simL, simR)     # Define the simulation box
simPath = dataPath + sim + '/' + simName + '/'

header = 'Column 1: Halo ID\n' + \
        'Column 2: Stellar mass (Msun)\n' + \
        'Column 3: Stellar half-mass radius (kpc)\n'

for snap in snapList:

    catalogue = load_snapshot(f'{simPath}SOAP/halo_properties_{snap:04d}.hdf5')
   
    halo_IDs = catalogue.input_halos.halo_catalogue_index.value

    Mstar = unyt.unyt_array(catalogue.bound_subhalo.stellar_mass.to_physical())
    Rstar = unyt.unyt_array(catalogue.bound_subhalo.half_mass_radius_stars.to_physical())

    SEL = (Mstar >= unyt.unyt_quantity(1e9, 'Msun')) * (Mstar <= unyt.unyt_quantity(1.05e9, 'Msun')) # Simple stellar mass selection. Replace this with 
    # your selection criteria.

    print(len(SEL[SEL]), 'galaxies selected in snapshot', snap)

    sample_file = np.vstack((halo_IDs, Mstar.to('Msun').value, Rstar.to('kpc').value)).T[SEL, :]

    np.savetxt(sampleFolder + 'sample_' + str(snap) + '.txt', sample_file, fmt = ['%d', '%.6e', '%.4f'], header = header)