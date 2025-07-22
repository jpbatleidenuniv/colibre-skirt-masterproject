"""
This script saves the SKIRT inputs for the galaxies in a given snapshot.
It reads the halo IDs and stellar masses from a sample file, loads the corresponding snapshot and SOAP catalogue,
and computes the necessary parameters for the SKIRT inputs.

Required is that not one single halo produces input files with just one particle"""

import os

import numpy as np
import tools.parse_sim as parse_sim
from swiftgalaxy.halo_catalogues import SOAP
from swiftgalaxy.iterator import SWIFTGalaxies, SWIFTGalaxy
from swiftsimio import load as load_snapshot
from swiftsimio.visualisation import generate_smoothing_lengths as gsl
from unyt.array import unyt_array, unyt_quantity

args, config = parse_sim.parse_sim()

starforming_gas_header, old_stars_header, dust = (
    config["Headers"]["starforming_gas"],
    config["Headers"]["old_stars"],
    config["Headers"]["dust"],
)
old_stars_tmin = unyt_quantity(
    config["SelectionCriteria"]["oldStarsAge"]["tmin"], "Myr"
)
preload_fields = config["Preload"]["preload_fields"]


skirt_runs_dir = config["SKIRTPaths"]["skirtRunsPath"]
skirt_inputs_dir = os.path.join(skirt_runs_dir, "inputs")
if not os.path.exists(skirt_inputs_dir):
    os.makedirs(skirt_inputs_dir, exist_ok=True)
else:
    print(f"Directory {skirt_inputs_dir} already exists, will overwrite existing files")


def compute_sml(sg: SWIFTGalaxy, Mstar: unyt_quantity) -> unyt_array:
    """Compute the smoothing lengths for the stars based on the total stellar mass of a galaxy."""
    if Mstar >= unyt_quantity(10 ** (8.5), "Msun"):
        # For large stellar masses, use a more accurate smoothing length calculation
        return (
            gsl(
                (sg.stars.coordinates + sg.centre) % sg.metadata.boxsize,
                sg.metadata.boxsize,
                kernel_gamma=1.0,
                neighbours=65,
                speedup_fac=2,
                dimension=3,
            )
            .to("pc")
            .to_physical()
        )
    else:
        # For smaller stellar masses, use the existing smoothing lengths
        return (
            sg.stars.smoothing_lengths.to("pc").to_physical() * 2.018932
        )  # Convert to physical units


def gas(sg):
    """For making the starforming gas particles, we use the TODDLERS SEDFamily with a normalized SFR. With the BPASS stellar template
    Requirements:
    - The averaged star formation rate in 10 Myr must be provided i.e. > 0. Msun/yr
    """
    # Star forming gas parameters
    x, y, z = sg.gas.coordinates.to("pc").value.T
    sml = sg.gas.smoothing_lengths.to("pc").to_physical() * 2.018932
    Z = sg.gas.metal_mass_fractions.to_physical()
    SFE = unyt_array(
        np.full(len(sml), 0.025), "dimensionless"
    )  # Star formation efficiency
    n_cloud = unyt_array(np.full(len(sml), 320.0), "1/cm**3")  # Cloud density
    SFR_averaged = (
        sg.gas.averaged_star_formation_rates[:, 1].to_physical().to("Msun/yr")
    )  # 10-Myr averaged SFRs

    star_forming_gas = np.transpose(
        [x, y, z, sml, Z, SFE, n_cloud, SFR_averaged]
    )  # Return the parameters in the required format
    # Dust paramaters (just the dust masses)
    M = sg.gas.masses.to("Msun").to_physical()
    grain_types = sg.gas.dust_mass_fractions.named_columns
    f_dust = unyt_array(
        [
            sg.gas.dust_mass_fractions.__getattribute__(type).to_physical()
            for type in grain_types
        ]
    )  # has shape (6, N)
    dust_masses = M * f_dust
    # print(f"Dust masses: {dust_masses} \n Shape: {dust_masses.shape}")
    dust = np.transpose(
        [x, y, z, sml] + [dust_masses[i] for i in range(len(grain_types))]
    )  # Return the dust parameters in the required format
    return np.atleast_2d(star_forming_gas), np.atleast_2d(dust)


def old_stars(sg: SWIFTGalaxy, Mstar):
    """For making the old stars (specified in SKIRTconfig.yaml), we use the BPASS SEDFamily using a Chabrier IMF where stars are not larger than 100 Msun.
    Requirements for BPASS Chab100 SEDFamily:
    - Stellar mass must be less than 100 Msun
    - Age must greater than 10 Myr
    - Need columns cooordinates, initial_mass, metallicity and age (6 columns in total)"""

    x, y, z = sg.stars.coordinates.to("pc").value.T
    sml = compute_sml(sg, Mstar)
    M = sg.stars.initial_masses.to("Msun").to_physical()
    Z = sg.stars.metal_mass_fractions.to_physical()
    age = sg.stars.ages.to("yr").to_physical()

    print(f"In total there are {age.shape} star particles.")

    # Now we filter for old stars
    age_mask = age >= old_stars_tmin
    if len(x) == 0 or not np.any(age_mask):
        return np.array([])  # If no old stars, return an empty array
    else:
        print(
            f"Found {np.sum(age_mask)} old stars with age >= {old_stars_tmin.value:.2f} Myr"
        )
        return np.atleast_2d(np.transpose([x, y, z, sml, M, Z, age]))[
            age_mask, :
        ]  # Return the parameters in the required format


def save_skirt_inputs(sg, halo_ID, Mstar):
    """Save the SKIRT inputs for a given halo ID and stellar mass."""
    # Create the filename for the SKIRT input file
    old_stars_params = old_stars(sg, Mstar)
    star_forming_gas_params, dust_params = gas(sg)

    #
    # If one of the paramater files contains no particles, we return an empty array
    # Save the parameters to the SKIRT input files

    np.savetxt(
        os.path.join(skirt_inputs_dir, f"ID{halo_ID}_old_stars.txt"),
        old_stars_params,
        header=old_stars_header,
        fmt="%.6e",
    )
    np.savetxt(
        os.path.join(skirt_inputs_dir, f"ID{halo_ID}_starforming_gas.txt"),
        star_forming_gas_params,
        header=starforming_gas_header,
        fmt="%.6e",
    )
    np.savetxt(
        os.path.join(skirt_inputs_dir, f"ID{halo_ID}_dust.txt"),
        dust_params,
        header=dust,
        fmt="%.6e",
    )
    print(
        f"Saved SKIRT inputs for halo ID {halo_ID} with stellar mass {Mstar:.2e} Msun"
    )


catalogue_file = config["ColibrePaths"]["haloCataloguePath"]
snap_file = config["ColibrePaths"]["virtualSnapshotPath"]
print(
    f"Loading snapshot {snap_file} and catalogue {catalogue_file} for snap {args.snap_nr}"
)

# Getting the information for the sample file
sample_file_path = os.path.join(skirt_runs_dir, f"sample_{args.snap_nr}.txt")
if not os.path.exists(sample_file_path):
    raise FileNotFoundError(
        f"Sample file {sample_file_path} does not exist. Please create it first."
    )
halo_IDs, Mstar = np.loadtxt(sample_file_path, usecols=[0, 2], unpack=True)
halo_IDs = halo_IDs.astype(int)  # Ensure halo IDs are integers
Mstar = unyt_array(
    Mstar, units="Msun", dtype=float
)  # Ensure Mstar is a unyt array with units

# Load the snapshot and SOAP catalogue
catalogue = load_snapshot(catalogue_file)
halo_IDs_all = catalogue.input_halos.halo_catalogue_index.value
halo_indices = np.where(np.isin(halo_IDs_all, halo_IDs))[0]
soap = SOAP(soap_file=catalogue_file, soap_index=halo_indices)

# We perform a final check to ensure that the loaded halos match in Halo IDs
if not np.array_equal(soap.input_halos.halo_catalogue_index.value, halo_IDs):
    raise ValueError(
        "The loaded SOAP halo indices do not match the expected halo IDs from the sample file."
    )
else:
    print(
        f"Loading {len(soap.input_halos.halo_catalogue_index.value)} halos from the SOAP catalogue for snapshot {args.snap_nr}"
    )

# Seems redundant, we have Mstar saved in the sample file.
# Now we still need Mstar for each halo, which we can get from the SOAP catalogue
# Mstar = unyt_array(soap.bound_subhalo.stellar_mass.to('Msun').to_physical())

# Making SWIFTGalaxy object
sgs = SWIFTGalaxies(
    snapshot_filename=snap_file, halo_catalogue=soap, preload=preload_fields
)

sgs.map(save_skirt_inputs, args=list(zip(halo_IDs, Mstar)))
print(f"Finished saving SKIRT inputs for snapshot {args.snap_nr}")


print(f"All SKIRT inputs saved in {skirt_inputs_dir}")
# End of the script
# Note: The script assumes that the sample file is already created and contains the halo IDs and stellar masses.
# The script also assumes that the SKIRT config file is correctly set up with the necessary paths and headers.
# The script will overwrite existing files in the skirt inputs directory.
