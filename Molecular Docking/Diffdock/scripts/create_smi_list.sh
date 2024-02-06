liganddir="/mnt/storage/valerie/Diffdock_docking/structures/ligand/Top100_Vina_results/all_sdf_files"

ligand_files=$liganddir/*

for l in $ligand_files
do
   while IFS= read -r line
   do
   	vars=( $line )
   	echo "'${vars[0]}'"
   done < "$l"
done
