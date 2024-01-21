# Molecular Docking with Glide 
The molecular docking was performed with Glide, a docking software developed by Schrödinger Inc.  
To access Glide and all other software developed by Schrödinger Inc., Maestro was used.
## Protein
As the structure for the receptor we chose the helicase with the ID 6ZSL.  
Important is that we only used the monomer version where only chain A is kept. The zinc ions were retained. 
The resulting pdb file was imported into Maestro and prepared using the included preparation tool. 
## Ligands
The top 100 ligands with the lowest affinity scores from the molecular docking with AutoDock Vina were used in the molecular docking with Glide. The files were imported into Maestro as sdf files.  
The ligands were also prepared using an included preparation tool.
## Receptor Grid
The same coordinates used to define the binding pocket in AutoDock Vina were also used to define the receptor grid for Glide:  
x = -14.466  
y = 12.57  
z = -73.214  
size = 30
## Results
|index|ZINC|ECBD|docking score|
|---|---|---|---|
|0|ZINC000096077632|EOS100380|-7.573|
|1|ZINC000008101127|nan|-7.273|
|2|ZINC000035880991|EOS100897|-6.786|
|3|ZINC000150588351|nan|-6.655|
|4|ZINC000072139230|EOS2443|-6.598|
|5|ZINC000014806405|EOS101245|-6.569|
|6|ZINC000085654587|EOS716|-6.554|
|7|ZINC000065362103|EOS2388|-6.552|
|8|ZINC000011679756|EOS100645|-6.549|
|9|ZINC000011679756|EOS100645|-6.549|
|10|ZINC000150338709|EOS100926|-6.517|
|11|ZINC000006745272|EOS101246|-6.41|
|12|ZINC000043208321|EOS100649|-6.215|
|13|ZINC000065232480|EOS1323|-6.092|
|14|nan|EOS100431|-6.089|
|15|nan|EOS101013|-6.054|
|16|nan|EOS100164|-5.988|
|17|ZINC000033266501|EOS815|-5.983|
|18|ZINC000150338819|nan|-5.972|
|19|ZINC000150338755|EOS101269|-5.937|