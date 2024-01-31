# 2023-2024_Team5-SU
Sorbonne University Team 5 

- Anna AUDIT
- Hortense BEAUSSANT
- Célia MESSAOUDI

Approach : 
Creating a basic pipeline to screen candidate ligands for a protein 

Pipeline :
- Generate 3D conformations from SMILES formula using RDKit (3 conformers for each candidate)
- Docking using AutoDock Vina
- Optimization of the results by doing molecular dynamics simulations (Making-it-rain) on the top scoring candidates

Files and folders :
- 'PDB_from_SMILES' : generates 3D conformers from SMILES formula of a ligand and saves them as PDB files
- 'conformations_ligands_pdb' : folder containing pdb files generated from SMILES formula (outputs of 'PDB_from_SMILES')
- 'Sorbonne5_flowchart' : png file summing up our pipeline
- 'docking_loop' : folder containing the code files to prepare ligands, protein and run and loop for all candidates
- 'find_top_scorers' : finds best scorers in folder of docking results
- 'top_docking_scores' : contains the list of the 29 best scoring compounds along with their score
- 'top_md_sim_scores' : contains table of ligands and their scores
- 'Sorbonne5_long_report' : long version of the report
- 'Sorbonne5_short_report' : chorter version of the report (2 pages)

![Team5_flowchart (2)](https://github.com/cu-bioinformatics/meet-eu-2023-projects/assets/148443412/a4ed48b7-4925-4b36-969e-0783498d00bf)

