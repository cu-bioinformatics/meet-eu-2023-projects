from Bio.PDB import PDBList
import subprocess
import pandas as pd
import os
from tqdm import tqdm 


def retrieve_nsp13(structure):
    """
    Download a NSP13 structure from the PDB.

    Parameters
    ----------
    structure : str
        The PDB identifier of the NSP13 structure to be downloaded.

    Returns
    -------
    str
        The content of the downloaded .cif file representing the structure.
    """

    pdbl = PDBList()
    pdbl.download_pdb_files([structure], pdir="NSP13")
    # .ent is another valid extension for pdb files. BioPython uses this extension
    os.makedirs("data/NSP13", exist_ok=True)
    fpath = f"data/NSP13/{structure}.cif"
    with open(fpath) as f:
        return f.read()
    

def generate_pocket():
    """
    Generate a pocket for the NSP13 structure using P2Rank.
    The generated pocket files will be saved in the 'pockets' folder.
    """
    # Path to the p2rank script (prank.sh)
    # should be in env variable
    p2rank_path = './p2rank/prank.sh'

    # Path to Git Bash for the Windows environment (to be replaced or improved if necessary)
    bash_path = 'c:/PROGRA~1/Git/bin/bash.exe'  # Make sure to set the correct path

    output_folder = 'data/pockets'
    os.makedirs("data/pockets", exist_ok=True)
    command = f'{bash_path} {p2rank_path} predict -f ./NSP13/6zsl.cif -o {output_folder}'
    subprocess.run(command, shell=True)


def get_pocket_coordinates(pocket_number):
    """
    Retrieve the coordinates of the specified pocket generated for the NSP13 structure.

    Parameters
    ----------
    pocket_number : int
        The number of the pocket to return coordinates.

    Returns
    -------
    tuple
        A tuple containing the X, Y, and Z coordinates of the pocket center.
    """
    pocket_path = 'data/pockets/6zsl.cif_predictions.csv'
    pockets_df = pd.read_csv(pocket_path)
    pockets_df.columns = [col.strip() for col in pockets_df.columns]
    best_pocket = pockets_df.query(f"rank == {pocket_number}") 
    (center_x, center_y, center_z) = best_pocket[["center_x", "center_y", "center_z"]].iloc[0]
    return (center_x, center_y, center_z)
    

def nsp13_pdb_to_pdbqt(structure_name):
    """
    Convert the NSP13 structure from PDB format to PDBQT format for docking.

    Parameters
    ----------
    structure_name : str
        The name of the NSP13 structure without file extension.
    """
    fpath_pdb = f"data/NSP13/{structure_name}.ent"
    fpath_pdbqt = f"data/NSP13/{structure_name}.pdbqt"
    
     # Convertis la liste de lignes en une seule chaîne
    with open(fpath_pdb) as f:
        pdb_content = f.read()

    pdb_content = '\n'.join(pdb_content)
    
    with open(fpath_pdb, 'w') as f:
        f.write(pdb_content)
    
    subprocess.run(['babel', fpath_pdb, '-O', fpath_pdbqt])