# 2023-2024_Team5-SU
Sorbonne University Team 5 
[this is not the final version]

Anna AUDIT

Hortense BEAUSSANT

Célia MESSAOUDI

Approach : 
Creating a basic pipeline to screen candidate ligands for a protein 

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
- 'top_docking_scores' : contains the list of the 29 best scoring compounds along with their score
- articles : reference articles

![image](https://github.com/cu-bioinformatics/meet-eu-2023-projects/assets/148443412/00850a1f-d2cf-4024-b862-04b416547f1a)
