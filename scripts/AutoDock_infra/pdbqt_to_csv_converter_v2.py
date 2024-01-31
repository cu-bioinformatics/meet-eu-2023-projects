import os, sys, argparse, time

def argsParser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--resultPath", required=True) # kde jsou XXX_output.pdbqt
    parser.add_argument("--outPath", required=True) # kde bude csv soubor
    args = vars(parser.parse_args())
    return args
args = argsParser()

RESULTS_PATH = args["resultPath"]
OUTPUT_PATH = args["outPath"] + "/" + time.strftime("%Y%m%d-%H%M%S") + ".csv" # csv soubot identifikovaný časem

path_generator = os.scandir(RESULTS_PATH)

if not os.path.exists(OUTPUT_PATH) or os.path.getsize(OUTPUT_PATH) == 0:
    with open(OUTPUT_PATH, "w") as outputfile:
        outputfile.write("proteinID,ligandID,pocketID,runID,affinity,rmsdLB,rmsdUB\n")
else: 
    print("error with output path")
    sys.exit(1)

for file_path in path_generator:
    if not file_path.name.startswith('.'):
        with open(file_path.path, "r") as file:
            # file name convention: receptor_ligand_pocket<n>_output.pdbqt
            fileIDs = file_path.name.split("_")
            proteinID = fileIDs[0]
            ligandID = fileIDs[1]
            pocketID = fileIDs[2][-1] # TODO just the pocket number

            with open(OUTPUT_PATH, "a") as outputfile:
                runID = 0
                for row in file:
                    if row[:19] == "REMARK VINA RESULT:" :
                        wantedIDs = row.split()
                        affinity = wantedIDs[3]
                        rmsdLB = wantedIDs[4]
                        rmsbUB = wantedIDs[5]
                        outputfile.write(str(proteinID) + ",")
                        outputfile.write(str(ligandID) + ",")
                        outputfile.write(str(pocketID) + ",")
                        outputfile.write(str(runID) + ",")
                        outputfile.write(str(affinity) + ",")
                        outputfile.write(str(rmsdLB) + ",")
                        outputfile.write(str(rmsbUB) + "\n")                    
                        runID += 1

