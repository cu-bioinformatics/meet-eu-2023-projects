#!/bin/bash

##############################################################################
dir_structures=/home/valerie/docking/structures
### receptor #######################################
receptor=/home/selina/docking/structures/protein/6zsl/monomer/6zsl_A.pdbqt
receptor_name=6zsl_monomer_a_2
### ligand #########################################
ligands=/home/valerie/docking/structures/ligand/ZINC/pdbqt/*.pdbqt
# ligands=/home/valerie/docking/structures/ligand/ECBD/pdbqt/*.pdbqt
### run version ####################################
version=ZINC
### results ########################################
dir_results=/home/selina/docking/results/$receptor_name
### config #########################################
# center
c_x=-14.466
c_y=12.57
c_z=-73.214
# size
s_x=30
s_y=30
s_z=30
# other
energy_range=3 #3
exhaustiveness=30 #8
num_modes=9 #9
cpu=28
seed=601399262

### create config txt file #################
# create txt file
config=conf.txt
touch $config
: > $config
# write txt file
# echo "ligand = ${ligand}" >> $config
echo "receptor = ${receptor}" >> $config
echo "center_x = ${c_x}" >> $config
echo "center_y = ${c_y}" >> $config
echo "center_z = ${c_z}" >> $config
echo "size_x = ${s_x}" >> $config
echo "size_y = ${s_y}" >> $config
echo "size_z = ${s_z}" >> $config
echo "energy_range = ${energy_range}" >> $config 
echo "exhaustiveness = ${exhaustiveness}" >> $config 
echo "num_modes = ${num_modes}" >> $config 
# echo "cpu = $cpu" >> $config
echo "seed = " $seed >> $config

    
version_file="version_${version}"

# create new version if config is different 
new=false
while [ $new = false ]
do 
    # if result directory of version already exists AND the new config file is different
    if [ -d $dir_results/$version_file ] && ! cmp -s $dir_results/$version_file/$config $config
    then 
        version=$((${version_file##*version_}+1))  # could lead to problems (DON'T CHANGE IT)
        echo $version
        version_file=${version_file%version_*}version_${version} 
    else
        # new results are created from different config file
        new=true
        dir_results=/home/selina/docking/results/$receptor_name/$version_file
        mkdir -p $dir_results/original_receptor
        
        # transfer relevant files to results
        cp $receptor $dir_results/original_receptor
        mv $config  $dir_results
        config=$dir_results/$config
        
    fi
done



##############################################################################
filter_file=size_filter.txt
# rm -f $dir_results/$filter_file
# molecule_size=ligand_size.txt

for ligand in $ligands
do
    ### ligand #################################
    ligand_name=${ligand##*/}
    ligand_name=${ligand_name%.*} # !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ### filter by size #########################
    scale=1
    correct_size=$(python3 GetMolecularSize.py $ligand $s_x $s_y $s_z $scale)
   if [ $correct_size = "True" ]
    then
    ### docking simulation ####################
        echo $ligand_name $correct_size 
        ./docking_simulation.sh $receptor_name $ligand $ligand_name $config $dir_results/ligand_results
    else
        echo $ligand_name $correct_size 
        touch $dir_results/$filter_file
        echo $ligand >> $dir_results/$filter_file
    fi 
done



