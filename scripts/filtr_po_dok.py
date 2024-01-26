"""filtr_po_dok-wer4.py"""

import argparse
import os
import re

# for plik in $(find *_G*/*.pdb); do python3 /media/mateusz/dysk_E/IPZ/filtr_po_dok.py ${plik} --output /media/mateusz/dysk_E/IPZ/wyniki_dokowania_baza/filtr_out_G.txt ; done
# for plik in $(find cabs_dock/*/*.pdb); do python3 /media/mateusz/dysk_E/IPZ/filtr_po_dok.py ${plik} --output /media/mateusz/dysk_E/IPZ/wyniki_dokowania_baza/filtr_out_cabs_dock.txt ; done
# for plik in $(find cabs_dock/*_G*/*.pdb); do python3 /media/mateusz/dysk_E/IPZ/filtr_po_dok.py ${plik} --output /media/mateusz/dysk_E/IPZ/wyniki_dokowania_baza/filtr_out_G.txt ; done
# for plik in $(find ./G*/*.pdb); do python3 /media/mateusz/dysk_E/IPZ/filtr_po_dok.py ${plik} --output /media/mateusz/dysk_E/IPZ/wyniki_dokowania_baza/filtr_out_G.txt ; done

parser = argparse.ArgumentParser(
        prog='Filtr po dokowaniu',
        description='znajdowanie sensownych dokowań',
        epilog='Text at the bottom of help')

parser.add_argument('input', nargs=1,help="plik pdb z samym peptydem lub z peptydem i helikazą ", default=None)
parser.add_argument('--bind',"-b", nargs=1,help="plik tekstowy z resztami wchodzącymi w skład kieszeni wiązania", default=["/media/mateusz/dysk_E/IPZ/binding_pockets_1.txt"])
parser.add_argument('--nsp13',"-n", nargs=1,help="plik pdb z helikazą nsp13", default=["/media/mateusz/dysk_E/IPZ/6zsl_A.pdb"])
parser.add_argument('--output',nargs=1, default=None)

args = parser.parse_args()

m=0

if os.path.isfile(args.input[0]):
	ins = args.input[0]
	f =open(ins,"r")
	line=f.readline().strip()
	if line.split()[0]=="ATOM":	
		m=1
	elif line[:14]=="REMARK Number:":
		m=2		
		if os.path.isfile(args.nsp13[0]):
			nsp13=args.nsp13[0]
		else:
			print("brak pliku nsp13")
	else:
		m=3
	f.close()

print(ins)
#print(m)

B={}

if args.bind:
	if os.path.isfile(args.bind[0]):
		with open(args.bind[0],"r") as b:
			for l,line in enumerate(b):
				B[l]=[]
				line=line.strip().split(",")
				for e in line:
					e=e.strip()
					if e:
						if m==1:
							if "relaxed" in ins:
								B[l].append([e[0],str(int(e[1:])+4)])
							else:
								B[l].append([e[0],str(int(e[1:])+2)])
						else:
							B[l].append([e[0],e[1:]])
	else:
		print(">> brak pliku z kieszenią wiązania\n")
if args.output is None:
	output="filtr_out.txt"
	f=open(output,"a")
	f.close()
elif os.path.isfile(args.output[0]):
	output = args.output[0]
elif os.path.exists(os.path.dirname(args.output[0])):
	output=args.output[0]
	f=open(output,"a")
	f.close()
else:
	output="filtr_out.txt"
	f=open(output,"a")
	f.close()
if os.stat(output).st_size==0:
	inp = open(output,"w")
	inp.write("""# 0 - AMP-binding pocket;
# 1 - 5'-binding pocket;
# 2 - kieszeń pomiędzy domeną wiążącą Zn, a stalk domain;
# 3 - reszta dokowań
""")

# list of missing residues in 6zsl_A
Missing=[95, 96, 97, 98, 99, 100, 101, 102, 186, 187, 188, 189, 190, 191, 192, 193,203,204,205,206,593,594,595,596,597,598,599,600,601]

S1=["A", "C","R","H","K","D","E", "S", "T", "N", "Q", "U", "G", "P", "V", "I", "L", "M", "F", "Y", "W"]
S3=["ALA", "CYS","ARG","HIS","LYS","ASP","GLU", "SER", "THR", "ASN", "GLN", "SEC","GLY", "PRO", "VAL","ILE","LEU", "MET", "PHE", "TYR", "TRP"]

Kod  =dict(zip(S1,S3))
ReKod=dict(zip(S3,S1))

###### !!!!!!! odległość, w której oczekujemy sasiadów
#### parametr odległość

o=7
		
def odl(a,b,d):
	w=0
#	print(f"a={a}, b={b}")
	for x,y in zip(a,b):
		w+=(x-y)**2
#	if w<=d**2:
#		print(f"\n>>w = {w}, d^2={d**2}\n")
	return w<=d**2
        
def fin(line, chain,pro="A"):
	line=line.split()
	if line[0]=="ATOM" and line[2]=="CA":
		if line[4]==chain:
			if (len(line)<12 and len(line)>9) or (m==2 and len(line)<9):
				spr=0
				for h in range(6,min(9,len(line))):
					if len(line[h])>8:
						spr=1
				if spr:
					if len(line[6])<=8:
						l6=line[6]
						l=line[7].split("-")
						if len(l)==2:
							l7=l[0]
							l8=f"-{l[1]}"
						else:
							l7=f"-{l[1]}"
							l8=f"-{l[2]}"
					else:
						l=line[6].split("-")
						if len(l)==2:
							l6=l[0]
							l7=f"-{l[1]}"
							l8=line[7]
						elif len(l)==3 and not l[0]:
							l6=l[0]
							l7=f"-{l[1]}"
							l8=f"-{l[2]}"
						elif len(l)==3 and l[0]:
							l6=f"-{l[1]}"
							l7=f"-{l[2]}"
							l8=line[7]	
						else:
							l6=f"-{l[1]}"
							l7=f"-{l[2]}"
							l8=f"-{l[3]}"	
					line[6]=l6
					line[7]=l7
					if len(line)==8:
						line.append(l8)
					else:
						line[8]=l8

			if chain==pro:
#				print("jestem tu")
#				print(line)
				Sym=ReKod[line[3]]
				Res=line[5]
				for i in range(len(B)):
					for r in B[i]:
						if r[0]==Sym and r[1]==Res:
							return [float(line[6]),float(line[7]),float(line[8]),i]
						elif r[1]==Res:
							print(float(line[6]),float(line[7]),float(line[8]),r[0], Sym)
			else:
#				print("jestem tu else")
				return [float(line[6]),float(line[7]),float(line[8]),-1]
				
#print(B)

if m==1:
####### tłumaczenie numeracji reszt do słownika B1
	print("Format z alphafold'a")
	w=[]
	wynik=0
	typ=[]
	pep=""
	pro=""
###### sprawdzamy, czy pierwszy, czy drugi łańcuch to peptyd.
	chains = open(ins,"r")
	line=1
	while not pep and line:
		line=chains.readline().strip().split()
		if len(line)==5:
			if int(line[4])<600:
				pep=line[3]
				
				if pep=="A":
					pro="B"
				else:
					pro="A"
#				print(pep)
#				print(pro)
	chains.close()
	with open(ins,"r") as chains:
		for line in chains:
			
			x =fin(line,pep,pro)
			if x:
				w.append(x)
			
	with open(ins,"r") as chains:
		for line in chains:
			x=fin(line,pro,pro)
#			print(x)
			if x:
				for r in w:
					if odl(r[:3],x[:3],o):
						wynik+=1
						typ.append(x[3])
					
#	print(w)
#	print(wynik)
#	print(typ)

		
			
	
elif m==2:
	w=[]
	wynik=0
	typ=[]
	with open(ins,"r") as chainB:
		for line in chainB:
			x =fin(line,"B")
			if x:
				w.append(x)
#	print("#"*40)
#	print(w)
	with open(nsp13,"r") as nsp:
		for line in nsp:
			x=fin(line,"A")
			if x:
				for r in w:
					if odl(r[:3],x[:3],o):
						wynik+=1
						typ.append(x[3])
						typ.append(x[3])
#	print(wynik)
#	print(typ)
elif m==3:
	w=[]
	wynik=0
	typ=[]
	with open(ins,"r") as chainB:
		for line in chainB:
			x =fin(line,"D")
			if x:
				w.append(x)
#	print("#"*40)
#	print(w)
	with open(ins,"r") as nsp:
		for line in nsp:
			x=fin(line,"A")
			if x:
#				print(line)
				for r in w:
					if odl(r[:3],x[:3],o):
						wynik+=1
						typ.append(x[3])
#	print(wynik)
#	print(typ)
if m:
	if wynik>=3:
		q=-1
		ost=-1
		for i in range(len(B)):
			if q<typ.count(i):
				q=typ.count(i)
				ost=i
		czy=True
#		print(ins)
		with open(output,"r") as f:
			for line in f:
				if ins in line and czy:
					czy=False
		if czy:
#######  elegancka gwiazdka w pliku przy 0 i 1
			if ost<2:
				v="\t*"
			else:
				v=""
####### ________________________________________
			f=open(output,"a")
			f.write(f"{ins}\t{ost}{v}\n")
			f.close()
	else:
		f=open(output,"a")
		f.write(f"{ins}\t{3}\n")
		f.close()

		

