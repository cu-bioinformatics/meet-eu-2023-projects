# Meet-EU 2023

Team members: Emilie Doan, Brenda Enriquez, Birsu Guzel and Isabelle Wu <br>
Supervisors : Juliana Silva Bernardes, Elodie Laine, Vaitea Opuu	

## The code

- molzip_adapted.py : contains the functions needed to predict with the Molzip method.
 This code has been adapted from https://github.com/daenuprobst/molzip.git

- features.ipynb : code allows you to check the relevance of the feature prediction.

- pdbbind.ipynb : code allows you to check the relevance of the dataset of LP_PDBBind.csv.

- main_molzip_result.ipynb : contains code that predicts the affinity of potential ligands for our problem.
  
- neural_network.ipynb : contains code that predicts separately,  affinity values and categorized kd values of the LP_PDBind.csv dataset.



## The data
The data needed to run our notebooks is in the data folder

- test_feat.csv : dataset from rdkit, needed to test our algorithm

- LP_PDBBind.csv : contains a ligand-pocket database that binds with a certain affinity which is in the 'value' column.
This dataset is from https://github.com/THGLab/LP-PDBBind/blob/master/dataset/LP_PDBBind.csv

- pilot_library.csv : contains more than 5000 potential ligands that could inhibit nsp13 that we want to predict the affinity.
This dataset is given by our school and their features are known.

- pocket.txt : the pockets we found for NSP 13

- data_PDB_Kd.csv : dataset from "LP_PDBBind.csv" that we process to have only some ligand-pocket


## Our goal

To find fixation classifier of ligands by using Molzip to reduce the size of our database.

Use different machine learning method with the idea of compress like Molzip.


## Molzip

Unconventional method based on the information theory. This method can predict different features relatively well based on the compression rate (Gzip). We accurately predict the solubility rate of different ligands and compared it to real solubility. We think that we can also use it to get other features as molecular weight.


## Neural Network

## Result

- result_molzip.csv : result from main_molzip_result.ipynb


