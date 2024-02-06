

liganddir="/mnt/storage/valerie/Diffdock_docking/structures/ligand/Top100_Vina_results/all_sdf_files"
ligand_files=$liganddir/*
outputdir="/mnt/storage/valerie/Diffdock_docking/structures/ligand/Top100_Vina_results/all_smi_files"


for l in $ligand_files
do
    ligand_id=${l%.*}           #remove .sdf
    ligand_id=${ligand_id##*/}  #remove everything else but ligand_id
    echo "Start converting $ligand_id into smi"
    obable -i $l -O $outputdir/$ligand_id.smi
    echo "End of converting $ligand_id into smi"
done
