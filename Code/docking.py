import os
import shutil
import json
from tqdm import tqdm
import subprocess
from process_nsp import get_pocket_coordinates
ROOT_PATH = 'c:/Users/alexi/OneDrive/Bureau/meetuStudent/Meet-U-2023-2024'
DOCKINGJOB_PATH = os.path.join(ROOT_PATH, "dockingjob")


def generate_config_file(center_x, center_y, center_z, fname_ligand):
    """
    Generate a configuration file for docking a ligand with NSP13.

    Parameters
    ----------
    center_x : float
        X-coordinate of the docking box center.
    center_y : float
        Y-coordinate of the docking box center.
    center_z : float
        Z-coordinate of the docking box center.
    fname_ligand : str
        The name of the ligand file for docking.
    """
    ligand_name = fname_ligand.split("_")[0]
    config = f"""receptor = 6zsl.pdbqt
ligand = {fname_ligand}

center_x = {center_x}
center_y = {center_y}
center_z = {center_z}

size_x = 20
size_y = 20
size_z = 20

out = docking_output/output_{ligand_name}.pdbqt
log = docking_output/logs_{ligand_name}.txt
num_modes = 10
energy_range = 4"""

    os.makedirs(DOCKINGJOB_PATH, exist_ok=True)

    file_path = f"{DOCKINGJOB_PATH}/config.txt"
    with open(file_path, "w") as file:
        file.write(config)


def run_docking(center_x, center_y, center_z, fname_ligand):
    """For the docking, the following files have to be in the folder dockingjob:
        -the protein, in format .pdbqt
        -the ligand, in format .pdbqt
        -the config file, in format.txt

    The output of the docking is saved at dockingjob/docking_output
    """
    # Remove the previous .pdbqt file
    previous_docking_ligand = [f for f in os.listdir(DOCKINGJOB_PATH) if f.endswith(".pdbqt") and "ZINC" in f]
    if previous_docking_ligand:
        fpath_ligand_to_remove = os.path.join(DOCKINGJOB_PATH, previous_docking_ligand[0])
        os.remove(fpath_ligand_to_remove)

    # Add the NSP13 to use for docking
    fpath_nsp13 = f"data/NSP13/6zsl.pdbqt"
    shutil.copy2(fpath_nsp13, DOCKINGJOB_PATH)
    
    # Add the ligand to use for docking
    fpath_ligand = f"data/ligands/zinc_db_pdbqt/{fname_ligand}"
    shutil.copy2(fpath_ligand, DOCKINGJOB_PATH)

    # Generate config file
    generate_config_file(center_x, center_y, center_z, fname_ligand)

    # Run docking
    os.makedirs(f"{DOCKINGJOB_PATH}/docking_output", exist_ok=True)
    vina_path = f"{DOCKINGJOB_PATH}/vina/vina.exe"
    subprocess.run([vina_path, "--config", "config.txt"], capture_output=False, text=True, cwd=DOCKINGJOB_PATH)


def run_docking_for_all_ligands(test=True, pocket_number=1):
    """
    Run docking for all prepared ligands with NSP13.

    Parameters
    ----------
    test : bool, optional
        Whether to run in test mode (default is True).
    pocket_number : int, optional
        The pocket number to use for docking (default is pocket 1).
    """
    ligands_folder = "data/ligands/zinc_db_pdbqt"
    fname_ligands = [f for f in os.listdir(ligands_folder) if f.endswith(".pdbqt")]

    # Always use the same pocket
    center_x, center_y, center_z = get_pocket_coordinates(pocket_number=pocket_number)
    for i, fname_ligand in tqdm(enumerate(fname_ligands), total=len(fname_ligands)):
        #print(f"Running docking for {fname_ligand.split('_')[0]}...")
        run_docking(center_x, center_y, center_z, fname_ligand)
        if test and i == 2:
            break


def retrieve_docking_results():
    """
    Retrieve docking results and save scores for each ligand in a JSON file.
    """
    all_scores = {}
    results_folder = f"{DOCKINGJOB_PATH}/docking_output"
    results_fname = [f for f in os.listdir(results_folder) if f.endswith(".pdbqt")]
    for fname in results_fname:
        ligand_name = fname.split("_")[1].split(".")[0]
        fpath_result = results_folder + "/" + fname
        with open(fpath_result) as f:
            result = f.readlines()[1]
            scores = [float(num) for num in result.split("RESULT:")[1].split() if float(num) != 0]
            all_scores[ligand_name] = max(scores)

    with open("output.json", "w") as f:
        json.dump(all_scores, f, indent=4)

    print("Score for each ligand saved at output.json")