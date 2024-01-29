#!/usr/bin/env dash

set -e

for i in lig/*.sdf; do

	smiles=$(obabel "$i" -osmi -n | sed -E "s:/:／:g;s:\\\:＼:g;s:[\t ]+$::g" )
	
	if [ -d "$smiles" ]; then
		echo "PDBQT file with SMILES $smiles already exists"
		exit 1
	fi
	
	mkdir "lig/$smiles"

	obabel "$i" -Otemp.sdf -h && mk_prepare_ligand.py -i "temp.sdf" -o "lig/$smiles/l.pdbqt" && rm temp.sdf
done
