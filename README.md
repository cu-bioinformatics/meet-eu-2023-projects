# 2023-2024_Team5-SU
Sorbonne University Team 5 
[this is no the final version]

Approach : 
Creating a basic pipeleine to screen candidate ligands for a protein 

Pipeline :
- Generate 3D conformations from SMILES formula using RDKit (3 conformers for each candidate)
- Docking using AutoDock Vina
- Optimization of the results by doing molecular dynamics simulations (Making-it-rain) on the top scoring candidates

Files and folders :
- 'PDB_from_SMILES' : generates 3D conformers from SMILES formula of a ligand and saves them as PDB files
- 'conformations_ligands_pdb' : folder containing pdb files generated from SMILES formula (outputs of 'PDB_from_SMILES')
- 'Team5_flowchart' : png file summing up our pipeline
- 'docking_loop' : folder containing the code files to prepare ligands, protein and run and loop for all candidates
- 'find_top_scorers' : finds best scorers in folder of docking results
- articles : reference articles
