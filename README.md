# NSP13-Ligand Docking Pipeline

This repository contains a Python script for performing molecular docking of ligands with the NSP13 protein structure. The docking pipeline involves several steps including downloading the NSP13 structure, generating pockets, preparing ligands, and running docking simulations.

## Prerequisites

Make sure to run pip install -r requirements.txt

Additionally, ensure that you have the following external tools available in your environment:

- Git Bash (for Windows users)
- OpenBabel
- AutoDock Vina
- RDKit
- P2Rank
- wget (for ligand download)

## Steps

### 1. Retrieving NSP13 Structure

The `retrieve_nsp13` function downloads the NSP13 structure from the Protein Data Bank (PDB) using its PDB identifier. The structure is saved in .cif format.

### 2. Generating Pockets

The `generate_pocket` function generates binding pockets on the NSP13 structure using P2Rank. The generated pockets are saved as CSV files.

### 3. Preparing Ligands

The `prepare_ligands` function downloads ligand structures from the ZINC database, adds hydrogens, and generates 3D conformations using RDKit. The prepared ligands are saved in .sdf format.

### 4. Converting Ligands to PDBQT

The `convert_ligands_to_pdbqt` function converts the prepared ligands from .sdf format to PDBQT format using OpenBabel. The converted ligands are saved in the 'ligands/zinc_db_pdbqt' folder.

### 5. Running Docking Simulations

The `run_docking_for_all_ligands` function performs docking simulations for all ligands with the NSP13 protein using AutoDock Vina. Docking is performed using the binding pockets generated earlier. The docking results are saved in 'dockingjob/docking_output' folder.

### 6. Retrieving Docking Results

The `retrieve_docking_results` function retrieves docking results, calculates scores for each ligand, and saves them in a JSON file named 'output.json'.

## Usage

1. Install the required dependencies.
2. Make sure external tools like Git Bash, OpenBabel, and AutoDock Vina are available in your environment.
3. Run the `docking_pipeline` function with appropriate parameters, such as the list of ZINC IDs for ligands to be docked.

## Notes

- Ensure proper file paths are set for external tools and input/output directories.
- Make sure to set the correct PDB identifier for the NSP13 structure.
- Adjust docking parameters in the configuration file generated for AutoDock Vina if needed.

