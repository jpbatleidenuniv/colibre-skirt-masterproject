import tools.parse_sim as parse_sim
import os
from multiprocessing import Pool
import numpy as np
from xml.etree import ElementTree as ET

args, config = parse_sim.parse_sim()

# Setting up the relevant paths
skirt_runs_dir = config['SKIRTPaths']['skirtRunsPath']
skirt_inputs_dir = os.path.join(skirt_runs_dir, 'inputs') # Where the SKIRT inputs (old stars, starforming gas, dust) are stored
skirt_outputs_dir = os.path.join(skirt_runs_dir, 'outputs') # Where the SKIRT outputs are stored
if not os.path.exists(skirt_outputs_dir):
    os.makedirs(skirt_outputs_dir, exist_ok=True)

Nprocesses = config['SimulationParameters']['Nprocesses']  # Number of processes to use for parallel processing
template_skifile = config['SKIRTPaths']['templateSkiFile']  # Path to the template .ski file for SKIRT

# Getting the information for the sample file
sample_file_path = os.path.join(skirt_runs_dir, f"sample_{args.snap_nr}.txt")
if not os.path.exists(sample_file_path):
    raise FileNotFoundError(
        f"Sample file {sample_file_path} does not exist. Please create it first."
    )
halo_IDs, max_dust_fraction = np.loadtxt(sample_file_path, usecols=[0, 3], unpack=True)
halo_IDs = halo_IDs.astype(int)  # Ensure halo IDs are integers

edits = config['SkiFileEdits']['edits']
template_skifile = config['SKIRTPaths']['templateSkiFile']
ski = ET.parse(template_skifile)

def make_edit(ski: ET.ElementTree, edits: list, halo_ID: str, max_dust_fraction: str):
    root = ski.getroot()
    vars = {"halo_ID": halo_ID, "max_dust_fraction": max_dust_fraction}

    for edit in edits:
        xpath = edit['xpath']

        # Try to find all matching elements
        elements = root.findall(xpath)
        if not elements:
            print(f"Warning: No elements found for {xpath}")
            continue

        for elem in elements:
            # Single attribute edit
            if "attribute" in edit:
                attr = edit['attribute']
                val = str(edit['value']).format(**vars)
                elem.set(attr, val)

            # Multiple attribute edit
            elif "attributes" in edit:
                for attr, val in edit['attributes'].items():
                    elem.set(attr, str(val).format(**vars))

    # Write to file once at the end
    skifile_path = os.path.join(skirt_runs_dir, f'ID{halo_ID}.ski')
    ski.write(os.path.join(skifile_path), encoding='utf-8', xml_declaration=True)
    print(f"Edited SKIRT file for halo {halo_ID} saved as ID{halo_ID}.ski")
    return skifile_path

def run_skirt(halo_ID: int, max_dust_fraction: float):
    """
    Function to run SKIRT for a given halo ID and max dust fraction.
    """
    skifile_path = make_edit(ski, edits, str(halo_ID), str(max_dust_fraction))
    
    # Construct the command to run SKIRT

    command = f"skirt -b -t 8 -i {skirt_inputs_dir} -o {skirt_outputs_dir} {skifile_path}"
    
    # Execute the command
    os.system(command)
    print(f"SKIRT run completed for halo {halo_ID}")

def main():
    """
    Main function to set up and run SKIRT simulations in parallel.
    """
    with Pool(Nprocesses) as pool:
        pool.starmap(run_skirt, zip(halo_IDs, max_dust_fraction))
        
if __name__ == "__main__":
    main()
