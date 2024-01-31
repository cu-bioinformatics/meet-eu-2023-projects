#!/usr/bin/python
# -*- coding: UTF-8 -*-

import argparse, os
import pandas

def argsParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--recPath", required=True, help="Path to folder with receptors")
    parser.add_argument("--confPath", required=True, help="Path to config files")
    parser.add_argument("--outPath", required=True, help="Output path")
    parser.add_argument("--ligPath", required=True, help="Path to folder with ligands")
    args = vars(parser.parse_args())
    return args
    
args = argsParser()
RECEPTORS_PATH = args["recPath"]
POCKET_NUMBERS=[2,3]

POCKETS_SIZES = {
            "x1": 20,
            "y1": 20,
            "z1": 20,
            "x2": 20,
            "y2": 20,
            "z2": 20
}
CONFIG_PATH = args["confPath"]

EXHAUSTIVENESS = 8
LIGANDS_PATH=args["ligPath"]
OUTPUT_PATH = args["outPath"]

path_gen_receptors = sorted((rec for rec in os.scandir(RECEPTORS_PATH) if rec.name.endswith(".pdbqt")),key=lambda receptors: receptors.name)
if len(path_gen_receptors) == 0:
    print("no receptors found")
    exit(1)


for i in range(len(path_gen_receptors)):
    file_path_receptor = path_gen_receptors[i]
    pockets = pandas.read_csv(RECEPTORS_PATH+"/structure.pdb_predictions.csv")
    path_gen_ligands = os.scandir(LIGANDS_PATH)
    keys = pockets.keys()  # weird spaces between keys
    for file_path_ligand in path_gen_ligands:
        for pocket in POCKET_NUMBERS:
            with open(CONFIG_PATH+r"/"+file_path_receptor.name[:-6]+"_"+file_path_ligand.name[:-6]+"_pocket"+str(pocket)+".txt","w") as file:
                file.write("receptor = "+file_path_receptor.path+"\n")
                file.write("ligand = "+file_path_ligand.path+"\n")
                file.write("center_x = "+str(pockets[keys[6]][pocket-1])+"\n")
                file.write("center_y = "+str(pockets[keys[7]][pocket-1])+"\n")
                file.write("center_z = "+str(pockets[keys[8]][pocket-1])+"\n")
                file.write("size_x = 20 \n")
                file.write("size_y = 20 \n")
                file.write("size_z = 20 \n")
                file.write("exhaustiveness = "+str(EXHAUSTIVENESS)+"\n")
                file.write("out = "+OUTPUT_PATH+r"/"+file_path_receptor.name[:-6]+"_"+file_path_ligand.name[:-6]+"_pocket"+str(pocket)+"_output.pdbqt\n")
                file.write("num_modes = 1") # vypíše jen ten nejlepší hit
