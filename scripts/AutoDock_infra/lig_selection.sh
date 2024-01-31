#!/bin/bash
# jako input dostane soubor, kde na kazdem radku bude jedna molekula
# vytvori novou slozku, jmeno definujeme
# z ligands_all vybere jen ty, ktere se jmenem shoduji

if [ $# -ne 2 ]; then
    echo "Usage [file with mol names] [new dir name]"
    exit 1
fi
MOL_FILE=$1
OUT_DIR="ligands/$2"

if [ ! -f $MOL_FILE ]; then
    echo "File $MOL_FILE does not exist"
    exit 1
fi

rm -r $OUT_DIR # pokud slozka existuje, vymazeme obsah rekurzivne
mkdir -p $OUT_DIR # tvorba slozky 

while IFS= read -r line; do
    filename=$(basename "$line")
    cp "ligands/ligands_all/$line.pdbqt" "${OUT_DIR}/${filename}.pdbqt" # soubor molekuly jenom prekopiruje
done < $MOL_FILE