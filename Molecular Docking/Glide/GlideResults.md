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
|index|ZINC|ECBD|docking score|smiles|
|---|---|---|---|---|
|0|ZINC000096077632|EOS100380|-7.573|CC[C@H](C)[C@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(=O)[C@@H](NC(=O)[C@H](CCCNC(=N)N)NC(=O)[C@@H](N)CC(=O)O)C(C)C)C(=O)N[C@@H](Cc1c[nH]cn1)C(=O)N1CCC[C@H]1C(=O)O|
|1|ZINC000008101127|nan|-7.273|CC1(C)C(C=C/C=C\C=C/C=C2/N(CCCCS(=O)(=O)O)c3ccc4ccccc4c3C2(C)C)=[N+](CCCCS(=O)(=O)O)c2ccc3ccccc3c21|
|2|ZINC000035880991|EOS100897|-6.786|Cc1cnc(Nc2ccc(F)cc2Cl)nc1-c1c[nH]c(C(=O)N[C@H](CO)c2cccc(Cl)c2)c1|
|3|ZINC000150588351|nan|-6.655|COC(=O)N[C@H](C(=O)N1CCC[C@H]1c1nc(-c2ccc3c(c2)O[C@@H](c2ccccc2)n2c-3cc3cc(-c4c[nH]c([C@@H]5CCCN5C(=O)[C@@H](NC(=O)OC)C(C)C)n4)ccc32)c[nH]1)C(C)C|
|4|ZINC000072139230|EOS2443|-6.598|Cc1ccc(NC(=O)N2CCN(c3cc(-n4ccnc4)ncn3)CC2)c(C)c1|
|5|ZINC000014806405|EOS101245|-6.569|NC(=O)c1ccc(-c2nc(-c3ccc4c(c3)OCCO4)c(-c3ccccn3)[nH]2)cc1|
|6|ZINC000085654587|EOS716|-6.554|Cc1ccc(-n2nc(N3CCN(C(=O)CCCN4C(=O)c5ccccc5C4=O)CC3)ccc2=O)cc1C|
|7|ZINC000065362103|EOS2388|-6.552|COc1cccc(NC(=O)N2CCN(c3cc(-n4nc(C)cc4C)ncn3)CC2)c1|
|8|ZINC000011679756|EOS100645|-6.549|CC1=NN(c2ccc(C)c(C)c2)C(=O)/C1=N\Nc1cccc(-c2cccc(C(=O)O)c2)c1O|
|9|ZINC000150338709|EOS100926|-6.517|CCc1cc(Nc2nccc(-c3c(-c4ccc(OC)c(C(=O)Nc5c(F)cccc5F)c4)nc4ccccn34)n2)c(OC)cc1N1CCC(N2CCN(S(C)(=O)=O)CC2)CC1|
|10|ZINC000006745272|EOS101246|-6.41|CNC(=O)c1cc(Oc2ccc(NC(=O)Nc3ccc(Cl)c(C(F)(F)F)c3)c(F)c2)ccn1|
|11|ZINC000043208321|EOS100649|-6.215|N#Cc1ccc(COC[C@H]2O[C@@H](n3c(NCc4ccc(Cl)c(Cl)c4)nc4c(N)ncnc43)[C@H](O)[C@@H]2O)cc1|
|12|ZINC000065232480|EOS1323|-6.092|O=C(Cn1cncc(-c2nc(-c3ccccc3)cs2)c1=O)Nc1ccc2c(c1)OCO2|
|13|nan|EOS100431|-6.089|O=C1O[C@]2(CC[C@@H](C(=O)Nc3ccn(-c4ccccc4F)n3)CC2)c2cnccc21|
|14|nan|EOS101013|-6.054|Cl.O=C(Nc1ccccc1-c1cn2c(CN3CCNCC3)csc2n1)c1cnc2ccccc2n1|
|15|nan|EOS100164|-5.988|C/C(=N\NC(=O)C(NC(=O)c1ccccc1)c1n[nH]c(=O)c2ccccc12)c1cccc(Br)c1|
|16|ZINC000033266501|EOS815|-5.983|O=C(NCc1cccnc1)C1CCN(c2nccn3nc(-c4ccccc4)cc23)CC1|
|17|ZINC000150338819|nan|-5.972|COC(=O)N[C@H](C(=O)N1CC2(CC2)C[C@H]1c1ncc(-c2ccc3c(c2)C(F)(F)c2cc(-c4ccc5nc([C@@H]6[C@H]7CC[C@H](C7)N6C(=O)[C@@H](NC(=O)OC)C(C)C)[nH]c5c4)ccc2-3)[nH]1)C(C)C|
|18|ZINC000150338755|EOS101269|-5.937|CC1(C)CCC(CN2CCN(c3ccc(C(=O)NS(=O)(=O)c4ccc(NCC5CCOCC5)c([N+](=O)[O-])c4)c(Oc4cnc5[nH]ccc5c4)c3)CC2)=C(c2ccc(Cl)cc2)C1|
|19|nan|EOS100913|-5.891|CC1=C(C(=O)Nc2cc3cn[nH]c3cc2F)C(c2ccc(C(F)(F)F)cc2)CC(=O)N1|