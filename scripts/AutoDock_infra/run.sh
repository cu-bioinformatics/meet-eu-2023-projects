#!/bin/bash
#PBS -q default@meta-pbs.metacentrum.cz
#PBS -l walltime=36:0:0
#PBS -l select=1:ncpus=4:mem=1gb:scratch_local=5gb

module add python36-modules-gcc
module add vina
module add parallel/20200322

DATADIR=/storage/praha1/home/larkus/meeteu
cd $DATADIR

# CONFIG FILE CREATION

REC_DIR=$r
LIG_DIR=$c
JOB_DIR=$j
# if [ $# -ne 3 ]; then # argument check
#     if [ $# -eq 1 ]; then
#         IFS=' ' read -r var1 var2 var3 <<< "$1"
#         REC_DIR=$var1
#         LIG_DIR=$var2
#         JOB_DIR=$var3
#     else
#         echo "Usage: $0 [receptor path] [ligands path] [job path]"
#         exit 1
#     fi
# else
#     REC_DIR=$1
#     LIG_DIR=$2
#     JOB_DIR=$3
# fi

CONF_DIR="$JOB_DIR/configs"
OUT_DIR="$JOB_DIR/vinaOut"
LOG_FILE="$JOB_DIR/log.txt"
VINA_LOG_FILE="$JOB_DIR/vina_log.txt"
SKIPPED_FILE="$JOB_DIR/skipped.txt"

mkdir -p "$JOB_DIR"
mkdir -p "$CONF_DIR"
> "$LOG_FILE"
> "$VINA_LOG_FILE"
> "$SKIPPED_FILE"

echo "$(date "+%d.%m.%Y-%H:%M:%S") INIT DONE" >> "$LOG_FILE" # LOG

python AutoDock_Config_Gen.py --recPath "$REC_DIR" --ligPath "$LIG_DIR" --confPath "$CONF_DIR" --outPath "$OUT_DIR" 
if [ $? -ne 0 ]; then
    echo "$(date "+%d.%m.%Y-%H:%M:%S") CONFIG GEN ERROR $?" >> "$LOG_FILE" # LOG
    exit 1
fi

echo "$(date "+%d.%m.%Y-%H:%M:%S") CONFIG GEN DONE" >> "$LOG_FILE" # LOG

# VINA EXECUTION
mkdir -p "$OUT_DIR"

for file in "$CONF_DIR"/* # pro každou konfiguraci
do 
    TMP_FILE="$OUT_DIR/tmp_output.pdbqt" # jedinný soubor se všemy výsledky runů
    echo "$(date "+%d.%m.%Y-%H:%M:%S") VINA run START - $file" >> "$LOG_FILE" # LOG    

    # zjistí jméno XXX_output.pdbqt z cfg souboru
    while IFS= read -r line; do
        if [[ $line == out* ]]; then
            # cesta k souboru výstupu Viny
            OUT_PDBQT="${line#*= }"
            break
        fi
    done < "$file"

    for run in {1..10}; do
        error_output=$(timeout 90s vina --config "$file" 2>&1 > "$VINA_LOG_FILE")
        if [ $? -eq 124 ]; then
            echo "$(date "+%d.%m.%Y-%H:%M:%S") VINA run $run TIMEOUT" >> "$LOG_FILE" # LOG
            echo "$file" >> "$SKIPPED_FILE"
            break
        elif [ -n "$error_output" ]; then
            error_msg=$(echo "$error_output" | grep -v '^\s*$')
            echo "$(date "+%d.%m.%Y-%H:%M:%S") VINA run $run ended with ERROR $error_msg" >> "$LOG_FILE" # LOG
            break
        else
            cat "$OUT_PDBQT" >> "$TMP_FILE" # appendujeme tmp soubor výstupem z Viny
        fi
    done
    
    cat "$TMP_FILE" >> "$OUT_PDBQT" # výsledky předáme zpátky do out souboru
    > $TMP_FILE # vymazání obsahu
    > $VINA_LOG_FILE # vymazání obsahu

    echo "$(date "+%d.%m.%Y-%H:%M:%S") VINA run END" >> "$LOG_FILE" # LOG
done
rm $TMP_FILE
echo "$(date "+%d.%m.%Y-%H:%M:%S") VINA DONE" >> "$LOG_FILE" # LOG

# CSV FILE CREATION
python pdbqt_to_csv_converter_v2.py --resultPath "$OUT_DIR" --outPath "$JOB_DIR"
if [ $? -ne 0 ]; then
    echo "$(date "+%d.%m.%Y-%H:%M:%S") CSV FILE ERROR $?" >> "$LOG_FILE" # LOG
    exit 1
fi
echo "$(date "+%d.%m.%Y-%H:%M:%S") CSV FILE DONE" >> "$LOG_FILE" # LOG

# CLEANUP
# smažeme všechno nepotřebné v JOB_DIR
rm -r "$CONF_DIR"
rm -r "$OUT_DIR"
rm "$VINA_LOG_FILE"
echo "$(date "+%d.%m.%Y-%H:%M:%S") CLEANUP DONE" >> "$LOG_FILE"