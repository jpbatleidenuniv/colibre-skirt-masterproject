"""
Run a set of SKIRT simulations for given halo indices.
Created by Andrea Gebek on 29.11.2024
"""

import numpy as np
import subprocess
from multiprocessing import Pool

# Set Paths

sampleFolder = '/Users/agebek/Downloads/' # Folder to the galaxy sample files
txtFilePath = '/Users/agebek/Downloads/' # Path to the COLIBRE particle .txt files
SKIRTinputFilePath = '/Users/agebek/Downloads/' # Path where the SKIRT input files will be stored

# Set list of snapshots to postprocess

snapList = [56, 123] # List of snapshots

Nprocesses = 3 # How many SKIRT simulations you want to run in parallel.
# Can also be set to one to run the SKIRT simulations serially.
# Note that each SKIRT simulation runs with 4 threads by default
# (this number can also be changed).

def preprocess(snapList):
    # Generate a list of SKIRT simulation names and run the necessary preprocessing steps

    skifilenames = []

    for snap in snapList:

        halo_IDs, Rstar = np.loadtxt(sampleFolder + '/sample_' + str(snap) + '.txt', unpack = True, usecols = [0, 2])
        halo_IDs = halo_IDs.astype(int)

        for idx, ID in enumerate(halo_IDs):

            skifilenames.append('snap' + str(snap) + '_ID' + str(ID))

            # Save SKIRT input files

            subprocess.run(['python', 'saveSKIRTinput.py', str(snap), str(ID), txtFilePath, SKIRTinputFilePath])

            # Edit ski files

            subprocess.run(['python', 'editSkiFile.py', str(snap), str(ID), str(Rstar[idx]), txtFilePath, SKIRTinputFilePath])

    return skifilenames



def runSKIRT(skifilename):

    # Run skirt

    subprocess.run(['skirt', '-t', '4', '-b', skifilename]) # Run SKIRT with 4 threads (that's apparently quite optimal)
    # The -b option reduces the verbosity of the log (but the saved log file still contains all logging information)

    return skifilename

def main():

    skifilenames = preprocess(snapList)

    with Pool(processes = Nprocesses) as pool:
        
        pool.map(runSKIRT, skifilenames)

if __name__=="__main__":

    main()