# colibre-skirt
This private, read-only repository contains our up-to-date framework for postprocessing COLIBRE simulations with the SKIRT 3D dust radiative transfer code. We expect that some settings/choices of the postprocessing pipeline are going to change over the coming weeks. To run SKIRT simulations, you will have to install [SKIRT](skirt.ugent.be) and follow the three steps outlined here.

## 1. SKIRT input files
SKIRT needs the information for the stars and dust particles in `.txt`files. At the moment, our pipeline works with eight `.txt`files per galaxy: One for the evolved stellar populations for star particles with ages above 30 Myr (idx_old_stars.txt). These star particles are modelled with the BPASS templates. Another `.txt`file is needed for the star-forming regions which are identified as the star particles with ages below 30 Myr (idx_young_stars.txt). These star particles are modelled with the TODDLERS templates, assuming subgrid dust attenuation which is not resolved in COLIBRE. The six remaining `.txt`files contain the information about the dust particles (one per species). The `saveSKIRTinput.py` script stores all these `.txt` files in the `SKIRTinputFiles` folder using the `swiftgalaxy` code from Kyle Oman (https://github.com/SWIFTSIM/swiftgalaxy) for a specific COLIBRE simulation and galaxy sample. Of course, you can also use your own preferred way to interface with COLIBRE data to generate these `.txt`files. Extra note: For some galaxies, most of the gas particles contain a negligible amount of dust. The `saveSKIRTinput.py` script stores only the gas particles with the biggest dust masses, until 99.9% of the total dust mass of the particles is reached. This can bring down the volume of the `.txt`files by a factor of a few, especially for bigger galaxies. Also, note that for galaxies with specific dust masses below 1e-5 the galaxy is assumed to be dust-free, such that no dust particle files are stored and the SKIRT simulation will be performed without dust (i.e. without medium).

## 2. The SKIRT configuration file
All settings of the SKIRT simulation are given in the `template_v3.0.ski`. For your own work you will probably need to adjust this file, e.g. to customize the output you want (i.e. spectra/images/datacubes). Feel free to reach out to us in case you have any questions on how to do this!

## 3. Running SKIRT for a set of galaxies
To iterate over the galaxies for which `.txt`input files were downloaded, run the `runSKIRT.py` script on a machine where SKIRT is installed. This will iterate over the galaxies with a specified number of SKIRT simulations in parallel, using the SKIRT parameters specified in `template_v3.0.ski`. The output of the SKIRT simulation (if you don't change the configuration file, this is simply a panchromatic SED for every galaxy) will be placed in the `SKIRToutputFiles` folder. Extra note: The runtime of SKIRT scales in many cases (especially as long as the number of photon packets is below approximately 1e8) with the number of cells in the spatial grid. We use an empirical grid refinement criterion based on the dust mass surface density of each simulated galaxy to have grids with appropriate resolutions. This leads to an average speedup of a SKIRT simulation by a factor of a few compared to a more conventional approach where the same grid refinement criterion is used for all galaxies.

## Installation instructions

`git clone https://andreagebek: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@github.com/andreagebek/colibre-skirt.git`
`cd colibre-skirt`
`mkdir SKIRTinputFiles SKIRToutputFiles`


## Contact
In case of questions or comments, please reach out to Andrea Gebek (andrea.gebek@ugent.be).