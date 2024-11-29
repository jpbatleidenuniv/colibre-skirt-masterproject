"""
Run a set of SKIRT simulations for given halo indices.
Created by Andrea Gebek on 29.11.2024
"""

import numpy as np
import subprocess
from multiprocessing import Pool
import os

# First, determine from the available SKIRT input files which halos to process
# This can of course easily replaced by e.g. a .txt file with a list of halo indices

halo_indices = 0


Nprocesses = 3 # How many SKIRT simulations you want to run in parallel.
# Can also be set to one to run the SKIRT simulations serially.

# Global settings

boxSize = 1e5 # In pc
Npp = int(10**7.5)
binTreeMaxLevel = 36
version = 'v3.0'

# Edit ski file

def runSKIRT(halo_index):

    skifilename = str(halo_index) + '_' + version + '.ski'
    skifilename_template = 'template_' + version + '.ski'

    subprocess.run(['cp', skifilename_template, skifilename]) # copy the skirt file for each galaxy

    # Check whether there is a dust .txt file. If not, run SKIRT without medium.

    if os.path.isfile('SKIRTinputFiles/' + str(halo_index) + '_FeSilicatesLarge.txt'):
        # Calculate max dust fraction based on particle data

        dust_x, dust_y, dust_z, dust_m_LargeFeSilicates = np.loadtxt('SKIRTinputFiles/' + str(halo_index) + '_FeSilicatesLarge.txt', unpack = True, usecols = [0, 1, 2, 4])
        dust_m_LargeMgSilicates = np.loadtxt('SKIRTinputFiles/' + str(halo_index) + '_MgSilicatesLarge.txt', usecols = 4)
        dust_m_LargeGraphite = np.loadtxt('SKIRTinputFiles/' + str(halo_index) + '_GraphiteLarge.txt', usecols = 4)
        dust_m_SmallFeSilicates = np.loadtxt('SKIRTinputFiles/' + str(halo_index) + '_FeSilicatesSmall.txt', usecols = 4)
        dust_m_SmallMgSilicates = np.loadtxt('SKIRTinputFiles/' + str(halo_index) + '_MgSilicatesSmall.txt', usecols = 4)
        dust_m_SmallGraphite = np.loadtxt('SKIRTinputFiles/' + str(halo_index) + '_GraphiteSmall.txt', usecols = 4)
        dust_m = dust_m_LargeFeSilicates + dust_m_LargeMgSilicates + dust_m_LargeGraphite + dust_m_SmallFeSilicates + dust_m_SmallMgSilicates + dust_m_SmallGraphite # In Msun

        dust_r = np.sqrt(dust_x**2 + dust_y**2 + dust_z**2) * 1e-3 # In kpc

        dustMasses_sorted = dust_m[np.argsort(dust_r)]

        idx_halfmass = np.max(np.argwhere((np.cumsum(dustMasses_sorted) / np.sum(dustMasses_sorted)) <= 0.5))

        dustHalfMassRadius = np.sort(dust_r)[idx_halfmass]

        dustHalfMass = (np.sum(dust_m) / 2.)

        SigmaDust = dustHalfMass / (np.pi * dustHalfMassRadius**2) # In solar masses / kpc*^2

        maxDustFraction = np.clip(10**(2.5 - 1.5 * np.log10(SigmaDust)), a_min = 10**(-6.5), a_max = 10**(-3.5))

        subprocess.run(['perl', '-pi', '-e', 's/maxLevel=\"0/maxLevel=\"' + str(binTreeMaxLevel) + '/g', skifilename])

        subprocess.run(['perl', '-pi', '-e', 's/FeSilicatesLarge.txt/' + str(halo_index) + '_FeSilicatesLarge.txt/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/MgSilicatesLarge.txt/' + str(halo_index) + '_MgSilicatesLarge.txt/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/GraphiteLarge.txt/' + str(halo_index) + '_GraphiteLarge.txt/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/FeSilicatesSmall.txt/' + str(halo_index) + '_FeSilicatesSmall.txt/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/MgSilicatesSmall.txt/' + str(halo_index) + '_MgSilicatesSmall.txt/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/GraphiteSmall.txt/' + str(halo_index) + '_GraphiteSmall.txt/g', skifilename])

        subprocess.run(['perl', '-pi', '-e', 's/minX=\"-0/minX=\"' + str(-boxSize / 2.) + '/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/maxX=\"0/maxX=\"' + str(boxSize / 2.) + '/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/minY=\"-0/minY=\"' + str(-boxSize / 2.) + '/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/maxY=\"0/maxY=\"' + str(boxSize / 2.) + '/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/minZ=\"-0/minZ=\"' + str(-boxSize / 2.) + '/g', skifilename])
        subprocess.run(['perl', '-pi', '-e', 's/maxZ=\"0/maxZ=\"' + str(boxSize / 2.) + '/g', skifilename])

        subprocess.run(['perl', '-pi', '-e', 's/maxDustFraction=\"0/maxDustFraction=\"' + str(maxDustFraction) + '/g', skifilename])

    else:
        # Change the skirt simulation to noMedium

        subprocess.run(['perl', '-pi', '-e', 's/simulationMode=\"DustEmission/simulationMode=\"NoMedium/g', skifilename])

        with open(skifilename, 'r+') as fp:
            # read and store all lines into list
            lines = fp.readlines()
            # move file pointer to the beginning of a file
            fp.seek(0)
            # truncate the file
            fp.truncate()

            # start writing lines
            # iterate line and line number
            for number, line in enumerate(lines):
                if number <= 40 or number >= 197:
                    fp.write(line)


    subprocess.run(['perl', '-pi', '-e', 's/numPackets=\"0/numPackets=\"' + str(Npp) + '/g', skifilename])

    subprocess.run(['perl', '-pi', '-e', 's#SKIRTinputFiles#SKIRTinputFiles/' + sim_name + '#g', skifilename])

    subprocess.run(['perl', '-pi', '-e', 's/old_stars.txt/' + str(halo_index) + '_old_stars.txt/g', skifilename])
    subprocess.run(['perl', '-pi', '-e', 's/young_stars.txt/' + str(halo_index) + '_young_stars.txt/g', skifilename])

    subprocess.run(['perl', '-pi', '-e', 's/radius=\"0/radius=\"' + str(boxSize / 2.) + '/g', skifilename])


    # Run skirt

    subprocess.run(['skirt', '-t', '16', '-b', skifilename]) # Run SKIRT with 16 threads (that's apparently quite optimal)

    # Remove unneeded SKIRT output and move SEDs to output folder

    subprocess.run(['rm', skifilename, skifilename[:-4] + '_parameters.xml', skifilename[:-4] + '_log.txt']) # Note that this also removes the log file, which you might not want
    subprocess.run(['mv', str(halo_index) + '_' + version + '_SED_sed.dat', 'SKIRToutputFiles/' + sim_name + '/' + str(halo_index) + '_SED.dat'])

    return skifilename

def main():

    with Pool(processes = Nprocesses) as pool:
        
        pool.map(runSKIRT, halo_indices)

if __name__=="__main__":

    main()