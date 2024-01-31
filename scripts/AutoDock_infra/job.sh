#!/bin/bash
# job.sh <name> <job_dir> <protein> <ligands> <nJobs>
# TODO chci tam mít i wall time atp nebo to bude u vseho stejne?

if [ $# -ne 5 ]; then
    echo "Usage $0 [job name] [work dir name] [protein name] [ligands dir] [n of jobs]"
    exit 1
fi


START=1
END=$5
JOB_NAME=$1
WORK_DIR_NAME=$2
PROT_NAME=$3
LIGS_PATH=$4
for (( i=$START; i<=$END; i++)); do
    echo "qsub -N "${JOB_NAME}_$i" -v r=proteins/$PROT_NAME,c=$LIGS_PATH/$i,j=$WORK_DIR_NAME/$PROT_NAME/$i run.sh"
    qsub -N "${JOB_NAME}_$i" -v r=proteins/$PROT_NAME,c=ligands/$LIGS_PATH/$i,j=$WORK_DIR_NAME/$PROT_NAME/$i run.sh
done