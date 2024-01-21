#!/bin/bash

receptor_name=$1 #string
ligand=$2 #path to file
ligand_name=$3 #string
config=$4 #path to file
dir_results=$5 #path to directory

echo $config

# naming stuff
# date=$(date '+%Y_%m_%d')
docking_of=${receptor_name}_${ligand_name} 
output_filename=${ligand_name}_conformations.pdbqt
log_filename=Docking_${receptor_name}_${ligand_name}.txt

mkdir -p $dir_results/$ligand_name
mkdir -p $dir_results/$ligand_name/pdbqt
mkdir -p $dir_results/$ligand_name/logfiles

# autodock vina 
vina --ligand $ligand \
--config $config \
--out $dir_results/$ligand_name/pdbqt/$output_filename \
--log $dir_results/$ligand_name/logfiles/$log_filename

python3 extract_table_from_logfiles.py \
$dir_results/$ligand_name/logfiles/$log_filename
