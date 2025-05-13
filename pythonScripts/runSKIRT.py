"""
Run a set of SKIRT simulations for given halo indices.
Created by Andrea Gebek on 29.11.2024
"""

import numpy as np
import subprocess
from multiprocessing import Pool
import yaml
import argparse
import os

parser = argparse.ArgumentParser(
    description="Run a set of SKIRT simulations for given halo indices."
)

# Set simName if needed for output files
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
    help="Simulation mode (default: Thermal_non_equilibrium).",
)

parser.add_argument(
        "--nproc",
        type=int,
        default=3,
        help="Number of SKIRT simulations you want to run in parallel. Note that each SKIRT simulation runs with 4 threads by default.",
)

args = parser.parse_args()

sim = 'L{:04.0f}N{:04.0f}'.format(args.BoxSize, args.NumParticles) # L0025N0752, L0100N1504, L0250N0752, L0250N1504, L0500N1504
simName = sim + '/' + args.mode # adds the simulation mode to the simulation name, e.g. L0100N1504/Thermal_non_equilibrium

args = parser.parse_args()

# Define filepaths from parameter file
dir_path = os.path.dirname(os.path.realpath(__file__))
with open(f'{dir_path}/../SKIRT_parameters.yml','r') as stream:
    params = yaml.safe_load(stream)

simPath = params['ColibreFilepaths']['simPath'].format(simName=simName)
sampleFolder = params['SkirtFilepaths']['sampleFolder'].format(sim=sim) # Folder to the galaxy sample files
txtFilePath = params['SkirtFilepaths']['storeParticlesPath'].format(sim=sim) # Path to the COLIBRE particle .txt files
SKIRTinputFilePath = params['SkirtFilepaths']['SKIRTinputFilePath'].format(sim=sim) # Path where the SKIRT input files will be stored
SKIRToutputFilePath = params['SkirtFilepaths']['SKIRToutputFilePath'].format(sim=sim) # Path where the SKIRT output files will be stored

# Set list of snapshots to postprocess

Nprocesses = args.nproc

def preprocess(snapList):
    # Generate a list of SKIRT simulation names and run the necessary preprocessing steps

    skifilenames = []

    for snap in snapList:

        halo_IDs = np.loadtxt(sampleFolder + '/sample_' + str(snap) + '.txt', unpack = True, usecols = 0)
        halo_IDs = halo_IDs.astype(int)

        for idx, ID in enumerate(halo_IDs):

            skifilenames.append( 'snap' + str(snap) + '_ID' + str(ID) )

            # Save SKIRT input files

            subprocess.run(['python', f'{dir_path}/saveSKIRTinput.py', str(snap), str(ID), txtFilePath, SKIRTinputFilePath])

            # Edit ski files

            subprocess.run(['python', f'{dir_path}/editSkiFile.py', str(snap), str(ID), txtFilePath, SKIRTinputFilePath, simPath])

    return skifilenames



def runSKIRT(skifilename):

    # Run skirt

    subprocess.run(['skirt', '-t', '16', '-b', skifilename]) # Run SKIRT with 4 threads (that's apparently quite optimal)
    # The -b option reduces the verbosity of the log (but the saved log file still contains all logging information)

    return skifilename

def postprocess(snapList):

    # Get the SKIRT output files and move them to the output folder

    for snap in snapList:

        halo_IDs = np.loadtxt(sampleFolder + '/sample_' + str(snap) + '.txt', unpack = True, usecols = 0).astype(int)

        for idx, ID in enumerate(halo_IDs):

            sim_name = 'snap' + str(snap) + '_ID' + str(ID)

            # subprocess.run(['rm', sim_name + '.ski']) # Remove the SKIRT input file

            subprocess.run(['mv', sim_name + '_parameters.xml', SKIRToutputFilePath + sim_name + '_parameters.xml'])
            subprocess.run(['mv', sim_name + '_log.txt', SKIRToutputFilePath + sim_name + '_log.txt'])

            if os.path.isfile(sim_name + '_conv_convergence.dat'): # Check if the file exists (it only does if there is a SKIRT medium)
                subprocess.run(['mv', sim_name + '_conv_convergence.dat', SKIRToutputFilePath + sim_name + '_conv_convergence.dat'])
                
            subprocess.run(['mv', sim_name + '_lum_luminosities.dat', SKIRToutputFilePath + sim_name + '_lum_luminosities.dat'])

            # We use a frame instrument which outputs a different file format
            subprocess.run(['mv', sim_name + '_image_total.fits', SKIRToutputFilePath + sim_name + '_image_total.fits'])
    


def main():

    skifilenames = preprocess(args.snaps)

    with Pool(processes = Nprocesses) as pool:
        
        pool.map(runSKIRT, skifilenames)

    postprocess(args.snaps)

if __name__=="__main__":

    main()