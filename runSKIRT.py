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

parser.add_argument(
    "Resolution",
    type=int,
    help="Particle mass resolution of the simulation in log10(M/Msun).",
)

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
    default="Thermal", # Thermal AGN feedback with non-equilibrium chemistry
    help="Simulation mode (default: Thermal).",
)


parser.add_argument(
        "--nproc",
        type=int,
        default=3,
        help="Number of SKIRT simulations you want to run in parallel. Note that each SKIRT simulation runs with 4 threads by default.",
)

args = parser.parse_args()

sim = 'L{:03.0f}_m{:01.0f}'.format(args.BoxSize, args.Resolution)
simName = sim + '/' + args.mode

args = parser.parse_args()

# Define filepaths from parameter file
dir_path = os.path.dirname(os.path.realpath(__file__))
with open(f'{dir_path}/SKIRT_parameters.yml','r') as stream:
    params = yaml.safe_load(stream)

simPath = params['InputFilepaths']['simPath'].format(simName=simName)
sampleFolder = params['OutputFilepaths']['sampleFolder'].format(simPath=simPath) # Folder to the galaxy sample files
txtFilePath = params['OutputFilepaths']['storeParticlesPath'].format(simPath=simPath) # Path to the COLIBRE particle .txt files
SKIRTinputFilePath = params['OutputFilepaths']['SKIRTinputFilePath'].format(simPath=simPath) # Path where the SKIRT input files will be stored

print(sampleFolder)

# Set list of snapshots to postprocess

Nprocesses = args.nproc

def preprocess(snapList):
    # Generate a list of SKIRT simulation names and run the necessary preprocessing steps

    skifilenames = []

    for snap in snapList:

        halo_IDs, Rstar = np.loadtxt(sampleFolder + '/sample_' + str(snap) + '.txt', unpack = True, usecols = [0, 2])
        halo_IDs = halo_IDs.astype(int)

        for idx, ID in enumerate(halo_IDs):

            skifilenames.append( params['OutputFilepaths']['skiFilepath'].format(simPath=simPath) + 'snap' + str(snap) + '_ID' + str(ID) )

            # Save SKIRT input files

            subprocess.run(['python', f'{dir_path}/saveSKIRTinput.py', str(snap), str(ID), txtFilePath, SKIRTinputFilePath])

            # Edit ski files

            subprocess.run(['python', f'{dir_path}/editSkiFile.py', str(snap), str(ID), str(Rstar[idx]), txtFilePath, SKIRTinputFilePath, simPath])

    return skifilenames



def runSKIRT(skifilename):

    # Run skirt

    subprocess.run(['skirt', '-t', '4', '-b', skifilename]) # Run SKIRT with 4 threads (that's apparently quite optimal)
    # The -b option reduces the verbosity of the log (but the saved log file still contains all logging information)

    return skifilename

def main():

    skifilenames = preprocess(args.snaps)

    with Pool(processes = Nprocesses) as pool:
        
        pool.map(runSKIRT, skifilenames)

if __name__=="__main__":

    main()