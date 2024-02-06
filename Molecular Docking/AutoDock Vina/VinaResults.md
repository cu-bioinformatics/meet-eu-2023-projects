# Molecular Docking with AutoDock Vina
## Database
Molecular Docking was performed with AutoDock Vina on ligands from two different databases.  
From the ZINC database only the fda approved ligands were considered and from the ECBD database only the Pilot library was included.  
Most of the ligands have assigned IDs from both databases.

## Receptor Protein
As the input structure for the receptor we chose the helicase with the ID 6ZSL.  
Important is that we only used the monomer version where only chain A is kept. We also removed the zinc ions due to AutoDock Vina which is not able to process those correctly. 
The resulting pdb was converted into pdbqt using AutoDockTools. 

## Binding Pocket
The binding pocket used for AutoDock Vina was the result of the Binding Pocket Prediction using Fpocket, FTMap and P2Rank.  
The result which were used for Vina:  
center_x = -14.466  
center_y = 12.57  
center_z = -73.214  
size_x = 30  
size_y = 30  
size_z = 30  

## Results
|rank|ZINC|ECBD|affinity|
|---|---|---|---|
|0|ZINC000150338709|EOS100926|-10.3|
|1|ZINC000006745272|EOS101246|-10.1|
|2|ZINC000043205799|EOS100157|-10.0|
|3|ZINC000011679756|EOS100645|-9.9|
|4|ZINC000058541155|EOS100491|-9.8|
|5|ZINC000095926668|EOS100499|-9.8|
|6|ZINC000006716957||-9.8|
|7|ZINC000015841814|EOS698|-9.8|
|8||EOS100286|-9.7|
|9|ZINC000003925087|EOS100030|-9.7|
|10|ZINC000004077976|EOS434|-9.7|
|11|ZINC000096170445|EOS100871|-9.6|
|12||EOS100700|-9.5|
|13|ZINC000058631343|EOS101275|-9.5|
|14|ZINC000095564436|EOS100548|-9.5|
|15||EOS100509|-9.4|
|16|ZINC000014954335|EOS101313|-9.3|
|17|ZINC000070466416|EOS101050|-9.3|
|18|ZINC000085654587|EOS716|-9.3|
|19|ZINC000095616580|EOS101166|-9.3|