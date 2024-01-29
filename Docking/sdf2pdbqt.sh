#!/usr/bin/env dash

set -e

if [ -z "$1" ]; then
	echo "$0 sdf-file"
	exit 1
fi

if [ ! -f "$1" ]; then
	echo "File $1 not found"
	exit 1
fi


smiles=$(obabel "$1" -osmi -n | sed -E "s:/:／:g;s:\\\:＼:g;s:[\t ]+$::g" )

if [ -f "$smiles.pdbqt" ]; then
	echo "PDBQT file with SMILES $smiles already exists"
	exit 1
fi

tmpfile="h_$1"
while [ -f "$tmpfile" ]; do
	tmpfile="_$tmpfile"
done

obabel "$1" -O"$tmpfile" -h && mk_prepare_ligand.py -i "$tmpfile" -o "$smiles.pdbqt" && rm "$tmpfile"
