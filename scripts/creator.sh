for i in K*.fasta; do cat $i > nsp13_$i; echo "" >> nsp13_$i; cat helicase_7re1.fasta >> nsp13_$i; done
for i in nsp13*.fasta; do
echo "#!/bin/bash" > run_${i%.*}.sh
echo "#PBS -N ${i%.*}" >> run_${i%.*}.sh
echo "#PBS -l nodes=1:ppn=8:gpus=2" >> run_${i%.*}.sh
echo "#PBS -l walltime=720:00:00" >> run_${i%.*}.sh
echo "#PBS -o user/af_peptides/stdout_${i%.*}" >> run_${i%.*}.sh
echo "#PBS -e user/af_peptides/stderr_${i%.*}" >> run_${i%.*}.sh

echo "python3 /home/apps/alphafold/2.2.0/run_singularity.py --fasta-paths user/af_peptides/seq/$i --model-preset multimer --data-dir user/alphafold_data --output-dir user/af_peptides/output" >> run_${i%.*}.sh
done
