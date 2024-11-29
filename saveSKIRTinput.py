"""
Script to create SKIRT input .txt files, based on Kyle Oman's
swiftgalaxy framework. As of November 2024, swiftgalaxy
needs to be installed from github (see https://swiftgalaxy.readthedocs.io/en/latest/getting_started/index.html).
Created by Andrea Gebek on 18.11.2024
"""

from swiftgalaxy.halo_catalogues import SOAP
import numpy as np
from swiftgalaxy.iterator import SWIFTGalaxies
from swiftsimio.objects import cosmo_array
from datetime import datetime
from swiftsimio.visualisation.smoothing_length.generate import generate_smoothing_lengths as gsl
from swiftsimio import load as load_snapshot
import unyt

startTime = datetime.now()

simulation = 'L050_m6/Thermal_non_equilibrium'

catalogue_file = '/cosma8/data/dp004/colibre/Runs/' + simulation + '/SOAP/halo_properties_0123.hdf5'
virtual_snapshot_file = '/cosma8/data/dp004/colibre/Runs/' + simulation + '/SOAP/colibre_with_SOAP_membership_0123.hdf5'


# Define the galaxy sample here

catalogue = load_snapshot(catalogue_file)

Mstar = unyt.unyt_array(catalogue.exclusive_sphere_100kpc.stellar_mass.to_physical()) # Convert the cosmo arrays to unyt arrays (without the "Physical" attribute).


SEL = (Mstar >= unyt.unyt_quantity(10**8.5, 'Msun')) # Selection based on stellar mass


halo_indices = np.where(SEL)[0]
halo_IDs = catalogue.input_halos.halo_catalogue_index.value[SEL]
Mstar = Mstar[SEL]
Mdust = unyt.unyt_array(catalogue.exclusive_sphere_100kpc.dust_large_grain_mass.to_physical()[SEL] + catalogue.exclusive_sphere_100kpc.dust_small_grain_mass.to_physical()[SEL])
specific_dust_masses = (Mdust / Mstar).to('dimensionless').value

soap = SOAP(catalogue_file, soap_index = halo_indices)

preload_fields = {'stars.coordinates', 'stars.metal_mass_fractions', 'stars.initial_masses', 'stars.ages',
                  'gas.coordinates', 'gas.smoothing_lengths', 'gas.masses',
                  'gas.dust_mass_fractions.GraphiteLarge', 'gas.dust_mass_fractions.MgSilicatesLarge', 'gas.dust_mass_fractions.FeSilicatesLarge',
                  'gas.dust_mass_fractions.GraphiteSmall', 'gas.dust_mass_fractions.MgSilicatesSmall', 'gas.dust_mass_fractions.FeSilicatesSmall'}

sgs = SWIFTGalaxies(virtual_snapshot_file, soap, preload = preload_fields)

def analysis(sg, halo_index, halo_ID, specific_dust_mass):
    # this function can also have additional args & kwargs, if needed
    # it should only access the pre-loaded data fields

    print('Saving txt files for halo index:', halo_index, 'with halo ID:', halo_ID)
    

    # Star particles
    #

    stars_x, stars_y, stars_z = sg.stars.coordinates.to('pc').to_physical().T
    # Recalculate stellar smoothing lengths, following COLIBRE tutorials
    stars_sml = gsl((sg.stars.coordinates + sg.centre) % sg.metadata.boxsize, sg.metadata.boxsize,
                    kernel_gamma = 2.0, neighbours = 65, speedup_fac = 2, dimension = 3).to('pc').to_physical()
    stars_Z = sg.stars.metal_mass_fractions.to_physical()
    stars_M = sg.stars.initial_masses.to('Msun').to_physical()
    stars_age = sg.stars.ages.to('yr').to_physical()
    stars_SFE = cosmo_array(np.full(len(stars_age), 0.025), 'dimensionless', comoving = False) # Star-formation efficienty, 2.5%
    stars_n_cl = cosmo_array(np.full(len(stars_age), 320.), 'cm**3', comoving = False) # Cloud density

    young_stars_mask = (stars_age.to('Myr').value <= 30.)
    old_stars_mask = (~young_stars_mask)

    old_stars_params = np.transpose([stars_x, stars_y, stars_z, stars_sml, stars_M, stars_Z, stars_age])[old_stars_mask, :]

    old_stars_header = 'column 0: x [pc]\n' + \
                'column 1: y [pc]\n' + \
                'column 2: z [pc]\n' + \
                'column 3: smoothing length [pc]\n' + \
                'column 4: initial stellar mass [Msun]\n' + \
                'column 5: metallicity [dimensionless]\n' + \
                'column 6: age [yr]\n'
    
    np.savetxt('SKIRTinputFiles/' + str(halo_index) + '_old_stars.txt', old_stars_params, fmt = '%.4e', header = old_stars_header)


    young_stars_params = np.transpose([stars_x, stars_y, stars_z, stars_sml, stars_age, stars_Z, stars_SFE, stars_n_cl, stars_M])[young_stars_mask, :]
    
    young_stars_header = 'column 0: x [pc]\n' + \
                'column 1: y [pc]\n' + \
                'column 2: z [pc]\n' + \
                'column 3: smoothing length [pc]\n' + \
                'column 4: age [yr]\n' + \
                'column 5: metallicity [dimensionless]\n' + \
                'column 6: star formation efficiency [dimensionless]\n' + \
                'column 7: cloud density [1/cm3]\n' + \
                'column 8: initial stellar mass [Msun]\n'
    
    np.savetxt('SKIRTinputFiles/' + str(halo_index) + '_young_stars.txt', young_stars_params, fmt = '%.4e', header = young_stars_header)


    # Gas/dust particles
    #

    if specific_dust_mass < 1e-5:
        
        print('Galaxy with ID', halo_ID, 'and index', halo_index, 'has a specific dust mass below 1e-5. Hence, no dust particles are stored for this galaxy.')
        
        return None
    
    gas_x, gas_y, gas_z = sg.gas.coordinates.to('pc').to_physical().T
    gas_sml = sg.gas.smoothing_lengths.to('pc').to_physical()
    gas_M = sg.gas.masses.to('Msun').to_physical()

    DustSpecies = sg.gas.dust_mass_fractions.named_columns
    gas_fDust = np.array([getattr(sg.gas.dust_mass_fractions, name) for name in DustSpecies]).T

    dust_M = (gas_fDust * np.atleast_1d(gas_M)[:, np.newaxis].repeat(6, axis = 1)).to('Msun').to_physical()

    # Remove dust particles with negligible dust mass contributions

    dust_M_allspecies = np.sum(dust_M, axis = 1)

    sortIndices = np.argsort(dust_M_allspecies.value)[::-1]

    cumulativeSum = np.cumsum(dust_M_allspecies[sortIndices])


    criterion = ((cumulativeSum / cumulativeSum[-1]).value >= 0.999) # Take the most massive dust particles until we reach 99.9% of the total dust mass

    criterion_index = np.min(np.argwhere(criterion))

    threshold_dust_M = dust_M_allspecies[sortIndices][criterion_index] # Threshold for the total dust mass of the particle

    dust_mask = (dust_M_allspecies.to('Msun').value >= threshold_dust_M)

    dust_params = np.transpose([np.atleast_1d(gas_x)[:, np.newaxis].repeat(6 ,axis = 1),
                                 np.atleast_1d(gas_y)[:, np.newaxis].repeat(6 ,axis = 1),
                                 np.atleast_1d(gas_z)[:, np.newaxis].repeat(6 ,axis = 1),
                                 np.atleast_1d(gas_sml)[:, np.newaxis].repeat(6, axis = 1),
                                 dust_M])[:, dust_mask, :]

    dust_header = 'column 0: x [pc]\n' + \
        'column 1: y [pc]\n' + \
        'column 2: z [pc]\n' + \
        'column 3: smoothing length [pc]\n' + \
        'column 4: dust mass [Msun]\n'
        
    for k, species in enumerate(DustSpecies):
        np.savetxt('SKIRTinputFiles/' + str(halo_index) + '_' + species + '.txt', dust_params[k], fmt = '%.4e', header = dust_header)

    return None

# map accepts arguments `args` and `kwargs`, passed through to function, if needed
sgs.map(analysis, args = list(zip(halo_indices, halo_IDs, specific_dust_masses)))


print('Elapsed time:', datetime.now() - startTime)
