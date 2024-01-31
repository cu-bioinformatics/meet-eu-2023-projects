
# ### Libraries



from Bio.PDB import PDBList
import subprocess
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import os
import shutil
import json


# ### 1. Download NSP13 structure



def retrieve_nsp13(structure):
    """Download a NSP13 structure from the PDB
    Returns
    -------
    structure : .cif file of the structure
    """
    pdbl = PDBList()
    pdbl.download_pdb_files([structure], pdir="NSP13", file_format="pdb")
    # .ent is another valid extension for pdb files. BioPython uses this extension
    fpath = f"./NSP13/pdb{structure}.ent"
    with open(fpath) as f:
        return f.read()


# ### 2. Convert protein to .pdbqt file


def nsp13_pdb_to_pdbqt(structure_name, pdb_content_lines):
    structure_pdb = f"./NSP13/{structure_name}.ent"
    structure_pdbqt = f"./NSP13/{structure_name}.pdbqt"
    
     # Convertis la liste de lignes en une seule chaîne
    pdb_content = '\n'.join(pdb_content_lines)
    
    with open(structure_pdb, 'w') as f:
        f.write(pdb_content)
    
    subprocess.run(['babel', structure_pdb, '-O', structure_pdbqt])
    
    return structure_pdbqt


# ### 3. Generate pocket with P2rank


def generate_pocket():

    # Path to the p2rank script (prank.sh)
    p2rank_path = './p2rank/prank.sh'

    # Path to Git Bash for the Windows environment (to be replaced or improved if necessary)
    bash_path = 'c:/PROGRA~1/Git/bin/bash.exe'  # Make sure to set the correct path

    # Folder to store pocket prediction results
    output_folder = './pockets'

    # Command to execute p2rank and predict pockets
    command = f'{bash_path} {p2rank_path} predict -f ./NSP13/6zsl.cif -o {output_folder}'

    # Execution of the command
    subprocess.run(command, shell=True)


# ### 4. Retrieve the coordinates of the pocket



def get_pocket_coordinates(pocket_number):
    pocket_path = './pockets/6zsl.cif_predictions.csv'
    pockets_df = pd.read_csv(pocket_path)
    pockets_df.columns = [col.strip() for col in pockets_df.columns]
    best_pocket = pockets_df.query(f"rank == {pocket_number}") 
    (center_x, center_y, center_z) = best_pocket[["center_x", "center_y", "center_z"]].iloc[0]
    return (center_x, center_y, center_z)


# 4.bis. Display several pockets of the 6ZSL comformation



csv = f'./pockets/6zsl.cif_predictions.csv'
r_csv = pd.read_csv(csv)

df = pd.DataFrame(r_csv)
df.head(2)


# ### 5. Download ligands from ZINC database



def download_ligands():
    downloader = f'./ligands-downloader/ligand_downloader.py'
    subprocess.run(['python', downloader], capture_output=True, text=True)


# ### 6. Preparation of the ligands with RDKIT



def prepare_ligands():
    csv_file = './ligands/ligands_downloader/csv_files/substances.csv'  
    df = pd.read_csv(csv_file,sep = ';')

    ligands_folder = './ligands/zinc_database/'
    os.makedirs(ligands_folder, exist_ok=True)

    for zinc_id in df['zinc_id']:
        sdf_file = f'{ligands_folder}/{zinc_id}.sdf'

        # Charger les molécules depuis le fichier SDF
        suppl = Chem.SDMolSupplier(sdf_file)
        ligands = [mol for mol in suppl if mol is not None]

        for mol in ligands:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)


        os.makedirs('./ligands/zinc_db_prepared/', exist_ok=True)
        modified_sdf_file = f'./ligands/zinc_db_prepared/{zinc_id}_prepared.sdf'
        w = Chem.SDWriter(modified_sdf_file)
        for mol in ligands:
            w.write(mol)
        w.close()


# ### 7. Convert ligands to .pdbqt with OpenBabel



def convert_ligands_to_pdbqt():
    ligands_folder = './ligands/zinc_db_prepared/'
    output_folder = './ligands/zinc_db_pdbqt/'
    os.makedirs(output_folder, exist_ok=True)

    for zinc_id in df['zinc_id']:
        sdf_file = f'{ligands_folder}/{zinc_id}_prepared.sdf'
        output_pdbqt_file = f'{output_folder}/{zinc_id}_prepared.pdbqt'

        # Open Babel for SDF to PDBQT
        subprocess.run(['babel', sdf_file, '-O', output_pdbqt_file, '-xh'])  # -xh ajoute des hydrogènes si nécessaire



# ### 8. Generate configuration file for docking with AutoDocks Vina



def generate_config_file(center_x, center_y, center_z, fname_ligand):
    """generate configuration files for a ligand
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

    os.makedirs("dockingjob", exist_ok=True)

    file_path = f"dockingjob/config.txt"
    with open(file_path, "w") as file:
        file.write(config)


# ### 9. Run molecular docking


def run_docking(center_x, center_y, center_z, fname_ligand):
    """For the docking, the following files have to be in the folder dockingjob:
        -the protein, in format .pdbqt
        -the ligand, in format .pdbqt
        -the config file, in format.txt

    The output of the docking is saved at dockingjob/docking_output
    """
    # Remove the previous .pdbqt file
    previous_docking_ligand = [f for f in os.listdir("dockingjob") if f.endswith(".pdbqt") and "ZINC" in f]
    if previous_docking_ligand:
        fpath_ligand_to_remove = os.path.join("dockingjob", previous_docking_ligand[0])
        os.remove(fpath_ligand_to_remove)
    
    # Add the ligand to use for docking
    fpath_ligand = f"ligands/zinc_db_pdbqt/{fname_ligand}"
    shutil.copy2(fpath_ligand, "dockingjob")

    # Generate config file
    generate_config_file(center_x, center_y, center_z, fname_ligand)

    # Run docking
    os.makedirs("dockingjob/docking_output", exist_ok=True)
    vina_path = "../Meet-U-2023-2024/dockingjob/vina/vina.exe"
    subprocess.run([vina_path, "--config", "config.txt"], capture_output=False, text=True, cwd="dockingjob")


def run_docking_for_all_ligands(test=True):
    ligands_folder = "ligands/zinc_db_pdbqt"
    fname_ligands = [f for f in os.listdir(ligands_folder) if f.endswith(".pdbqt")]

    # Always use the same pocket
    center_x, center_y, center_z = get_pocket_coordinates(pocket_number=1)
    for i, fname_ligand in enumerate(fname_ligands):
        print(f"Running docking for {fname_ligand.split('_')[0]}...")
        run_docking(center_x, center_y, center_z, fname_ligand)
        if test and i == 2:
            break


# ### 10. Retrieve docking results and score


def retrieve_docking_results():
    all_scores = {}
    results_folder = "dockingjob/docking_output"
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


# ### 11. Create an end-to-end docking pipeline


def docking_pipeline(structure="6zsl"):
    nsp13 = retrieve_nsp13(structure)
    nsp13_pdb_to_pdbqt(structure, nsp13)
    generate_pocket()

    download_ligands()
    prepare_ligands()
    convert_ligands_to_pdbqt()

    run_docking_for_all_ligands()
    retrieve_docking_results()

if __name__ == "__main__":
    docking_pipeline(structure="6zsl")
