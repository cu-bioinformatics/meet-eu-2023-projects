"""
filtr_po_dok v2
"""

import argparse
import os
import re

# for plik in $(find */*.pdb); do python3 path_to_script ${plik} ; done
# for plik in $(find cabs_dock/*/*.pdb); do python3 path_to_script ${plik} --output dir_to_output_text_file ; done


parser = argparse.ArgumentParser(
    prog='Filtr po dokowaniu',
    description='znajdowanie sensownych dokowań',
    epilog='Text at the bottom of help')

parser.add_argument('input', nargs=1, help="plik pdb z samym peptydem lub z peptydem i helikazą ", default=None)
parser.add_argument('--bind', "-b", nargs=1, help="plik tekstowy z resztami wchodzącymi w skład kieszeni wiązania",
                    default=["path_to_text_file_with_binding_pocket_data.txt"])
parser.add_argument('--nsp13', "-n", nargs=1, help="plik pdb z helikazą nsp13", default=["path_to_protein_model.pdb"])
parser.add_argument('--output', nargs=1, default=None)

args = parser.parse_args()

m = 0

if os.path.isfile(args.input[0]):
    ins = args.input[0]
    f = open(ins, "r")
    line = f.readline().strip()
    if line.split()[0] == "ATOM":
        m = 1
    elif line[:14] == "REMARK Number:":
        m = 2
        if os.path.isfile(args.nsp13[0]):
            nsp13 = args.nsp13[0]
        else:
            print("brak pliku nsp13")
    else:
        m = 3
    f.close()

print(ins)

B = {}

if args.bind:
    bind_file = args.bind[0]
    if os.path.isfile(bind_file):
        B = [[] for _ in range(sum(1 for _ in open(bind_file)))]
        with open(bind_file, "r") as b:
            for l, line in enumerate(b):
                line = [e.strip() for e in line.strip().split(",") if e]
                if m == 1:
                    B[l] = [[e[0], str(int(e[1:]) + 2)] for e in line]
                else:
                    B[l] = [[e[0], e[1:]] for e in line]
    else:
        print(">> brak pliku z kieszenią wiązania\n")

if args.output is None:
    output = "filtr_out.txt"
else:
    output = args.output[0] if os.path.isfile(args.output[0]) else args.output[0] if os.path.exists(
        os.path.dirname(args.output[0])) else "filtr_out.txt"

with open(output, "a") as f:
    if os.stat(output).st_size == 0:
        f.write("""# 0 - AMP-binding pocket;
# 1 - 5'-binding pocket;
# 2 - kieszeń pomiędzy domeną wiążącą Zn, a stalk domain;
# 3 - reszta dokowań
""")

# list of missing residues in 6zsl_A
Missing = [95, 96, 97, 98, 99, 100, 101, 102, 186, 187, 188, 189, 190, 191, 192, 193, 203, 204, 205, 206, 593, 594, 595, 596, 597, 598, 599,
           600, 601]

S1 = ["A", "C", "R", "H", "K", "D", "E", "S", "T", "N", "Q", "U", "G", "P", "V", "I", "L", "M", "F", "Y", "W"]
S3 = ["ALA", "CYS", "ARG", "HIS", "LYS", "ASP", "GLU", "SER", "THR", "ASN", "GLN", "SEC", "GLY", "PRO", "VAL", "ILE", "LEU", "MET", "PHE",
      "TYR", "TRP"]

Kod = dict(zip(S1, S3))
ReKod = dict(zip(S3, S1))

o = 7


def odl(a, b, d):
    w = 0

    for x, y in zip(a, b):
        w += (x - y) ** 2

    return w <= d ** 2


def fin(line, chain, pro="A"):
	line = line.split()

	if line[0] == "ATOM" and line[2] == "CA" and line[4] == chain:
		if 9 < len(line) < 12:
			spr = 0
			for h in range(6, 9):
				if len(line[h]) > 8:
					spr = 1

			if spr:
				if len(line[6]) <= 8:
					l6 = line[6]
					l = line[7].split("-")

					if len(l) == 2:
						l7, l8 = l[0], f"-{l[1]}"
					else:
						l7, l8 = f"-{l[1]}", f"-{l[2]}"
				else:
					l = line[6].split("-")

					if len(l) == 2:
						l6, l7, l8 = l[0], f"-{l[1]}", line[7]
					elif len(l) == 3 and not l[0]:
						l6, l7, l8 = l[0], f"-{l[1]}", f"-{l[2]}"
					elif len(l) == 3 and l[0]:
						l6, l7, l8 = f"-{l[1]}", f"-{l[2]}", line[7]
					else:
						l6, l7, l8 = f"-{l[1]}", f"-{l[2]}", f"-{l[3]}"

				line[6], line[7], line[8] = l6, l7, l8

			if chain == pro:
				Sym, Res = ReKod[line[3]], line[5]
				for i, b_chain in enumerate(B):
					for r in b_chain:
						if r[0] == Sym and r[1] == Res:
							return [float(line[6]), float(line[7]), float(line[8]), i]
						elif r[1] == Res:
							print(float(line[6]), float(line[7]), float(line[8]), r[0], Sym)
			else:
				return [float(line[6]), float(line[7]), float(line[8]), -1]




def process_chain_file(file_path, pep, pro, output_file, distance_threshold):
	w = []
	wynik = 0
	typ = []

	with open(file_path, "r") as chains:
		for line in chains:
			x = fin(line, pep, pro)
			if x:
				w.append(x)

	with open(file_path, "r") as chains:
		for line in chains:
			x = fin(line, pro, pro)
			if x:
				for r in w:
					if odl(r[:3], x[:3], distance_threshold):
						wynik += 1
						typ.append(x[3])

	if wynik:
		if wynik >= 3:
			q = max(typ.count(i) for i in range(len(B)))
			ost = typ.index(q) if q in typ else -1

			if ost != -1:
				with open(output_file, "r") as f:
					if not any(file_path in line for line in f):
						v = "\t*" if ost < 2 else ""
						with open(output_file, "a") as f:
							f.write(f"{file_path}\t{ost}{v}\n")
		else:
			with open(output_file, "a") as f:
				f.write(f"{file_path}\t{3}\n")

# Main part
if m == 1:
	print("Format z alphafold'a")

	pep, pro = "", ""
	with open(ins, "r") as chains:
		for line in chains:
			line = line.strip().split()
			if len(line) == 5 and int(line[4]) < 600:
				pep = line[3]
				pro = "B" if pep == "A" else "A"
				break

	process_chain_file(ins, pep, pro, output, o)

elif m == 2:
	process_chain_file(ins, "B", "A", output, o)

elif m == 3:
	process_chain_file(ins, "D", "A", output, o)
