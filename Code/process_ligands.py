import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem
import os


def download_ligands(zinc_id_list):
    """
    Download ligand structures from the ZINC database.

    Parameters
    ----------
    zinc_id_list : list (or any iterable)
        A list of ZINC IDs corresponding to the ligands to be downloaded.
    """
    # Creation of a folder for the ligands
    ligands_folder = 'data/ligands/zinc_database'
    os.makedirs(ligands_folder, exist_ok=True)

    # Downloading of the SDF Files from the ZINC databse
    for zinc_id in zinc_id_list:
        sdf_url = f'https://zinc.docking.org/substances/{zinc_id}.sdf'
        command = f'wget {sdf_url} -O {ligands_folder}/{zinc_id}.sdf'
        os.system(command)


def prepare_ligands():
    """
    Prepare ligands by adding hydrogens and generating 3D conformations.
    The prepared ligands will be saved in the 'ligands/zinc_db_prepared' folder.
    """
    zinc_db_folder = "data/ligands/zinc_databse"
    zinc_fnames = [f for f in os.listdir(zinc_db_folder) if f.endswith(".sdf")]

    for fname in zinc_fnames:
        ligand_name = fname.split(".")[0]
        sdf_file = os.path.join(zinc_db_folder, fname)

        # Load molecule
        suppl = Chem.SDMolSupplier(sdf_file)
        ligands = [mol for mol in suppl if mol is not None]

        # Add hydrogens and build 3d conformation
        for mol in ligands:
            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol)

        # Saved prepared ligands
        os.makedirs('data/ligands/zinc_db_prepared/', exist_ok=True)
        modified_sdf_file = f'data/ligands/zinc_db_prepared/{ligand_name}_prepared.sdf'
        w = Chem.SDWriter(modified_sdf_file)
        for mol in ligands:
            w.write(mol)
        w.close()


def convert_ligands_to_pdbqt():
    """
    Convert prepared ligands from SDF format to PDBQT format using OpenBabel.
    The converted ligands will be saved in the 'ligands/zinc_db_pdbqt' folder.
    """
    output_folder = 'data/ligands/zinc_db_pdbqt/'
    os.makedirs(output_folder, exist_ok=True)
    zinc_prepared_folder = "data/ligands/zinc_db_prepared"
    zinc_fnames = [f for f in os.listdir(zinc_prepared_folder) if f.endswith(".sdf")]
    ligand_name = fname.split(".")[0]

    for fname in zinc_fnames:
        sdf_file = os.path.join(zinc_prepared_folder, zinc_fnames)
        output_pdbqt_file = f'{output_folder}/{ligand_name}_prepared.pdbqt'

        # Convert from sdf to pdbqt with OpenBabel
        subprocess.run(['babel', sdf_file, '-O', output_pdbqt_file, '-xh'])  # -xh ajoute des hydrogènes si nécessaire