echo "Start Diffdock docking process"

smi_and_id_file="/mnt/storage/valerie/Diffdock_docking/structures/ligand/Top100_Vina_results/smi_and_id.txt"
proteinpath="/mnt/storage/valerie/Diffdock_docking/structures/protein/6zsl_A.pdb"
outputdir="/mnt/storage/valerie/Diffdock_docking/results/Top100_vina/diffdock_output"


while IFS="" read -r line
do
	vars=( $line )
   	ligand_smi=${vars[0]}
   	ligand_id=${vars[1]}
   	for l in $line
   	do
   	echo "Start docking with $ligand_id"
   	python -m inference --protein_path $proteinpath --ligand $ligand_smi --out_dir $outputdir --complex_name $ligand_id --save_visualisation --inference_steps 20 --samples_per_complex 40 --batch_size 10 --actual_steps 18 --no_final_step_noise 
   	echo "Done with $ligand_id"
   	done
done < $smi_and_id_file
