import subprocess
import os
import sys
import csv
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen
from rdkit.Chem import AllChem
from rdkit.Chem import MolFromPDBFile
from rdkit.Chem import MolToPDBFile
from rdkit.Chem import MolToPDBBlock
from typing import IO
from meeko import MoleculePreparation
from meeko import PDBQTWriterLegacy

def usage():
    print("python3 main.py [Path to VinaGPU]")
    exit(1)


if len(sys.argv) != 2:
    usage()

autodock_gpu_path = sys.argv[1]

res_file_print_outputs = "res.txt"
really_nice = "really_nice.txt"
def prepare_structures():
    for i in range(10,11):
        pdb_file = f"resources/structures/structure_{i}_cleaned.pdb"
        prepare_protein(i, pdb_file)



def prepare_ligand(smiles: str):
    """
    processes ligand and returns it in pdbqt format; needed for screening
    """
    ligand_mol = Chem.MolFromSmiles(smiles)
    ligand_mol = AllChem.AddHs(ligand_mol)
    AllChem.EmbedMolecule(ligand_mol)
    preparator = MoleculePreparation()
    mol_setups = preparator.prepare(ligand_mol)
    for setup in mol_setups:
        ligand_pdbqt, is_ok, error_msg = PDBQTWriterLegacy.write_string(setup)
        if is_ok:
            return ligand_pdbqt
    return error_msg




def prepare_protein(id: int, pdb_file):
    """
    processes protein structure for screening and saves it in pdbqt format
    """
    protein_mol = MolFromPDBFile(pdb_file, sanitize=False)
    print(protein_mol)

    protein_mol = Chem.RemoveHs(protein_mol)
    protein_mol = Chem.DeleteSubstructs(protein_mol, Chem.MolFromSmarts('[OH2]'))

    protein_mol = AllChem.AddHs(protein_mol)
    AllChem.EmbedMolecule(protein_mol)
    AllChem.ComputeGasteigerCharges(protein_mol)
    '''
    AllChem.MergeNonPolarHs(protein_mol)
    '''''
    s = MolToPDBBlock(protein_mol)
    with open(f"resources/structures/structure_{id}.pdbqt", "w") as f:
        f.write(s)

def create_list_of_drugs(csv_file_path: str):
    data_dict = {}
    with open(csv_file_path, 'r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        for row in csv_reader:
            # Assuming the 'index' column is the id
            data_dict[row['chembl_id']] = row['canonical_smiles']
    return data_dict


def generate_docking_command(receptor_file, flex_file, ligand_file, center_x, center_y, center_z, size_x, size_y, size_z, output_file, num_modes, gpu_id=0):
    vina_command = [
            autodock_gpu_path,
            "--receptor", receptor_file,
            "--ligand", ligand_file.strip(),
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--out", output_file,
            "--num_modes", str(num_modes),
    ]

    return " ".join(vina_command)


def perform_docking(receptor_file, ligand_file, output_file, gpu_id=0):
    flex_file = "flex"
    center_x = -17
    center_y = -5
    center_z = 1
    size_x = 15
    size_y = 20
    size_z = 20
    vina_command = generate_docking_command(receptor_file, flex_file, ligand_file, center_x, center_y, center_z, size_x, size_y, size_z, output_file, num_modes=3, gpu_id=0)
    try:
        subprocess.run(vina_command, shell=True)
        with open(res_file_print_outputs, "a") as f:
            f.write(f"Docking for {ligand_file} completed successfully.\n")
        return True
    except subprocess.CalledProcessError as e:
        with open(res_file_print_outputs, "a") as f:
            f.write(f"Error during docking for {ligand_file}: {e}\n")
        return False


def do_screening(id: str, smiles: str):
    """
    returns boolean; if screening succeeded
    """
    # prepare ligand
    try:
        ligand_pdbqt = prepare_ligand(smiles)
        with open(f"{id}.pdbqt", "w") as f:
            f.write(ligand_pdbqt)
        # screening
        ligand_file = f"{os.path.abspath(f'{id}.pdbqt')}\n"
        output_file = "output.pdbqt"
        for i in range(10,11):
            protein_file ="ex.pdbqt"  # f"resources/structures/structure_{i}.pdbqt"
            output_file = "output.pdbqt"
            t_f = perform_docking(protein_file, ligand_file, output_file, gpu_id=0)
            try:
                successful_docking = analyze_docking_results(id, output_file, "src/docking_results.tsv")
                if successful_docking:
                    break
            except:
                with open("fatal.txt", "a") as f:
                    f.write(f"{id},{smiles}\n")
                os.remove(output_file)
                return False
        # remove temporary files
        os.remove(f"{id}.pdbqt")
        return successful_docking
    except:
        return 0

def analyze_docking_results(id: str, output_file: str, tsv_file):
    energy, rmsd1, rmsd2 = extract_energy_rmsd(output_file, id)
    #if float(rmsd1) < 6:
    append_result_to_tsv(id, energy, rmsd1, rmsd2, tsv_file)
    return float(rmsd1) < 3


def extract_energy_rmsd(output_file: str, id):
    with open(output_file, 'r') as f:
        lines = f.readlines()

    energy = None
    rmsd1 = None
    rmsd2 = None

    for line in lines:
        if "REMARK VINA RESULT:" in line:
            energy, rmsd1, rmsd2 = line.split(":")[1].strip().split()
            if float(energy) < -5.0:
                with open(really_nice, "a") as f:
                    f.write(f"{id} : {energy}\n")
            break

    return float(energy), float(rmsd1), float(rmsd2)

def append_result_to_tsv(id, energy, rmsd1, rmsd2, tsv_file):
    df = pd.DataFrame({'Ligand id': [id], 'Free energy of binding (kcal/mol)': [energy], 'RMSD (Angstroms)': [rmsd1], 'RMSD2': [rmsd2]})
    df.to_csv(tsv_file, mode='a', header=False, index=False, sep='\t')


def screening(path_to_files_csv):
    drugs = create_list_of_drugs(path_to_files_csv)

    for id, smiles in drugs.items():
        do_screening(id, smiles)


if __name__ == '__main__':
    #prepare_structures()
    screening("potential_drugs.csv")