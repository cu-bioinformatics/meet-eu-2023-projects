import pandas as pd
import os
from openbabel import openbabel

def convert_sdf_to_smiles(sdf_path):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats("sdf", "smi")

    mol = openbabel.OBMol()
    if obConversion.ReadFile(mol, sdf_path):
        return obConversion.WriteString(mol).strip()
    else:
        return None

def process_ligands(csv_path, sdf_dir, output_csv):
    # Read CSV file
    df = pd.read_csv(csv_path)
    smiles_data = []

    for ligand_id in df['ligandID']:
        sdf_file = os.path.join(sdf_dir, f"{ligand_id}.sdf")

        if os.path.exists(sdf_file):
            smiles = convert_sdf_to_smiles(sdf_file)
            if smiles:
                smiles_data.append({'LigandID': ligand_id, 'SMILES': smiles})
            else:
                print(f"Failed to convert {ligand_id}")
        else:
            print(f"SDF file for {ligand_id} not found.")

    # Create new DataFrame and save to CSV
    smiles_df = pd.DataFrame(smiles_data)
    smiles_df.to_csv(output_csv, index=False)
    print(f"Output written to {output_csv}")

# Replace 'input.csv' and 'sdf_directory' with your actual file path and directory
process_ligands('/Users/faflik/Projects/MoleculeScreen/Data/results/repurposing/7RDY/updated_repurposing_result.csv', '/Users/faflik/Projects/MoleculeScreen/Conversion/SdfToPdbtq/struc_files', 'output_smiles.csv')
