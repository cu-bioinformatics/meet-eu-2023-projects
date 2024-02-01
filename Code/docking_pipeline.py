import pandas as pd
from process_nsp import retrieve_nsp13, nsp13_pdb_to_pdbqt, generate_pocket
from process_ligands import download_ligands, prepare_ligands, convert_ligands_to_pdbqt
from docking import run_docking_for_all_ligands, retrieve_docking_results


def docking_pipeline(zinc_id_list, structure="6zsl", pocket_number=1, test=True):
    """
    Perform the entire docking pipeline for NSP13-ligand docking.

    Parameters
    ----------
    zinc_id_list : list
        A list of ZINC IDs corresponding to the ligands to be docked.
    structure : str, optional
        The PDB identifier of the NSP13 structure (default is "6zsl").
    pocket_number : int, optional
        The number of the pocket to use for docking (default is 1).
    test : bool, optional
        Whether to run in test mode (default is True).
    """
    retrieve_nsp13(structure)
    nsp13_pdb_to_pdbqt(structure)
    generate_pocket()

    download_ligands(zinc_id_list)
    prepare_ligands()
    convert_ligands_to_pdbqt()

    run_docking_for_all_ligands(pocket_number=pocket_number, test=test)
    retrieve_docking_results()

if __name__ == "__main__":
    zinc_id_list = pd.read_csv("data/substances.csv", sep=";")["zinc_id"]
    docking_pipeline(zinc_id_list, test=True)