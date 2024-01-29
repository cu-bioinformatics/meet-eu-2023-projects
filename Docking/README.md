# Docking

Place your SMILES.txt file in this folder


## Preparing files
First of create SDF files (if you alredy have sdf files in folder you can copy them to lig directory)

     python create_files.py
## Create folders
With files created use 

    ./makeligdirs.sh

to create correct directories structure for docking script
## Running Vina
Then run Vina

    ./dock.sh
 To generate summary of docking in /lig folder run
 

    ./generate_tsvs.sh


