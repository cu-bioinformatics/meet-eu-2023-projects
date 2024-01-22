from rcsbsearchapi.search import StructSimilarityQuery, SeqMotifQuery
import requests
import json
from Bio.PDB import PDBList
import os
import re

pdb_list = PDBList()

peptydy=["SVFSGYRK", "CLANGMIMY", "RQQPGKGPRY", "ALEATCKSL", "LDAQSAPLRV",
	"FLRQNEVL", "IDALNENK", "DIEQLRSQL", "IQKVAGTW", "TLPFHSVIY", "VLEKYKDVI", "LVGPTIWR",
	"RFPRPVY", "DFVCRAFQY", "YLYQGPIVL", "HMWPGDIK", "ALKIPISKIY", "LTFQHNF", "ALTEGLPEL", "VLMKDVQEI",
	"KLTCNLTR", "WVRQAPGKALE", "SYFSGLDPYL", "IIVTQTMK", "LWGYTDPFL", "ILGFHANE", "AIKIPLL",			# , "QIPTVNNL" - nie chce wczytać .pdb
	"ALQDIIGSL", "ILFYVKEF", "ITGFYPEE", "GLKAVFPLL", "EIIESPLF", "GYFYPIQI", "LSITENGEFK", "GWLEPLL",
	"SLPYPFI", "YVHPFHL", "GLLPGLMVY", "GLDIQK", "HLSPNDPIF", "MIQLDLI", "SFVSEVPEL", "GLPFPPEL", "LHTPLPL"]
	
path_pdb="dir_to_pdb"
path_from_pdb="dir_to_pdb"
path_to_save_fasta="dir_to_candidates"
path_to_K_fasta="dir_to_fastas"

"""
funkcja fasta_z_pdb() pobiera sekwencje fasta dla zadanych struktur .pdb peptydów n, zapisanych 
pod nazwa {n}_motif.pdb w katalogu pf 
"""

def fasta_z_pdb(n,fp):
	q4 = StructSimilarityQuery(structure_search_type="file_upload",
                           file_path=f"{fp}/{n}_motif.pdb",  # specify local model file path
                           file_format="pdb")
	list(q4())


	last=""
	with open(f"{path_to_save_fasta}/{n}.fasta","w") as f:
		for i in range(len(q4("polymer_entity").iquery())):
			pdb_id=q4("polymer_entity").iquery()[i]
			last=pdb_id.split("_")
			nr=int(last[1])
			last=last[0]
			req = requests.get(f"https://www.rcsb.org/fasta/entry/{last}/display").text
			l=req.strip().split("\n")
			if len(l[2*nr-1])<=2*len(n):
				pep=f"{l[2*nr-1]}"
				Xs=0
				Xm=0
				Xe=0
				if pep[0]=="X":
					pep1=re.sub(r"^X","",pep)
					Xs=len(pep)-len(pep1)
					pep=pep1
				if pep[-1]=="X":
					pep1=re.sub(r"X$","",pep)
					Xe=len(pep)-len(pep1)
					pep=pep1
				if "X" in pep:
					pep1=pep.replace("X","")
					Xm=len(pep)-len(pep1)
					pep=pep1
				f.write(f"{l[2*nr-2]}\tXs={Xs}\tXm={Xm}\tXe={Xe}\n{pep}\n")
						
						
#fasta_z_pdb(n,path_pdb)					
						



"""
Znajduje PDB_id i pozycje motywu dla reprezentatywnych białek ze znalezionej grupy
Chimera oparła się uruchomieniu skryptu.
"""

def pdb_z_motif(n,path_from_pdb):
#	w=[]
	my_query = {
		"query": {
			"type": "terminal",
			"service": "seqmotif",
			"parameters": {
				"value": f"{n}",
      			"pattern_type": "simple",
				"sequence_type": "protein"
			}
		},
		"return_type": "polymer_entity",
		"request_options": {
			"results_verbosity": "verbose",
			"group_by_return_type": "representatives",
			"group_by": {
				"aggregation_method": "sequence_identity",
				"ranking_criteria_type": {
					"sort_by": "rcsb_entry_info.resolution_combined",
					"direction": "asc"
				},
				"similarity_cutoff": 100
			},
			"paginate": {
				"start": 0,
				"rows": 25
			},
			"results_content_type": [
				"experimental"
			],
			"sort": [
				{
					"sort_by": "score",
					"direction": "desc"
				}
			],
			"scoring_strategy": "combined"
		}
	}
	
	my_query = json.dumps(my_query)
	data = requests.get(f"https://search.rcsb.org/rcsbsearch/v2/query?json={my_query}")
	if data.status_code==200:
		data=data.json()
		pdb_set=data['result_set']
		for pdb in pdb_set:
			ids=pdb["identifier"]
			dane=pdb['services'][0]['nodes'][0]['match_context']
			for s in dane:
#				w=(ids.split("_")[1],s['start'], s['end'])
#				filename=pdb_list.retrieve_pdb_file(ids.split("_")[0],pdir=path_from_pdb, file_format="pdb")

#				print(filename)
				print(f"{ids.replace('_','.')} start={s['start']}, end = {s['end']}")
			

#pdb_z_motif(n,path_from_pdb)
k=1

for p in peptydy:
###### znajdowanie plików reprezentatywnych z pdb:

#	pdb_z_motif(p,path_from_pdb)

###### znajdowanie fasta z pdb na podstawie stworzonych w chimerze plików pdb.
	if os.path.isfile(f"{path_pdb}/{p}_motif.pdb"):
		print(p)
		fasta_z_pdb(p,path_pdb)
		if os.stat(f"{path_to_save_fasta}/{p}.fasta").st_size:
			with open(f"{path_to_save_fasta}/{p}.fasta","r") as f:
				for line in f:
					if line[0]==">":
						g=open(f"{path_to_K_fasta}/K{k}.fasta","w")
						line=line.strip().split("\t")
						g.write(f">K{k}\t{line[1]}\t{line[2]}\t{line[3]}\n")
					else:
						g.write(line)
						g.close()
						k+=1
