echo -e "Start preparation of ligand ... \n"

read -p "Enter the directory of the file, where all your sdf files are stored. It probably starts with:  " dir

#Checking wether directory exists
until [ $dir ] && [ -d $dir ]
do
    read -p "Sorry, directory doesn't exist! Try again:  " dir
done

echo "directory exists"

#load all sdf files in dir
ligand_files=$dir/*

read -p "Enter the directory of the file, where you want to store your pdbqt files:  " outputdir

#Checking wether directory exists
until [ $outputdir ] && [ -d $outputdir ]
do
    read -p "Sorry, directory doesn't exist! Try again:  " ouptdir
done

echo "directory exists"

# repeat prepare_ligand.py for all sdf files in directory
for l in $ligand_files
do
    ligand_id=${l%.*}           #remove .sdf
    ligand_id=${ligand_id##*/}  #remove everything else but ligand_id
    # add -h to make hydrogens explicit for ZINC, but happens beforehand with ECBD
    #--gen3d: generate 3D structure
    # explicit hydrogens and 3d structure are needed for prepare_ligand.py from ADFR_suite
    # you might need to mkdir for files with exlicit hydrogens beforhand
    obabel -isdf "$l" -osdf -O "$dir/with_hydrogens/$ligand_id.with_hydrogens.sdf" -h --gen3d
    echo "obable done with $ligand_id"
    #call prepare_ligand from ADFR suite 
    mk_prepare_ligand.py -i "$dir/with_hydrogens/$ligand_id.with_hydrogens.sdf" -o "$outputdir/$ligand_id.pdbqt"
    echo "meeko is done with $ligand_id"
done

                                                              39,10         Bot


