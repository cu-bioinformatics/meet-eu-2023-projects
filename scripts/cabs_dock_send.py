import requests
import json
import argparse
import os
from Bio import SeqIO



parser = argparse.ArgumentParser(
        prog="Skrypt uruchamiający Cabs-dock'a\njob's name: {6ZSLA}_{file name}")

parser.add_argument('input', nargs=1,help="plik fasta z samym peptydem", default=None)
args = parser.parse_args()

url = 'https://biocomp.chem.uw.edu.pl/CABSdock/REST/add_job/'
# files = {'file': open('your_PDB_file.pdb')} #or use PDB code in var data

if os.path.isfile(args.input[0]):
	with open(args.input[0],"r") as handle:
		for record in SeqIO.parse(handle, "fasta"):
			seq=record.seq
	name=str(args.input[0]).split(".")[0]

data = {
	"receptor_pdb_code": "6ZSL:A", #or use PDB file in var files
	"ligand_seq": str(seq),
	"email": "mc430767@students.mimuw.edu.pl",
	"show": True,
	"project_name":f"6ZSLA_{name}",
	
}

#response = requests.post(url, files=files, data=data) #request with file
response = requests.post(url, data=data) # request without file
print(response.text)
