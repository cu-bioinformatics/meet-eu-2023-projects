################################################ 
import pandas as pd
from collections import Counter
import os
from pathlib import Path
import re
import shutil
################################################ 
# input 
################################################ 
top_n=""
receptor="6zsl_monomer_a"
print("This script takes the vina results and finds the top scores.\nPLEASE give the number of top scorers you want.")
while top_n=="":
    top_n = input("Number of top scorers:")
    try:
        top_n=int(top_n)
    except:
        print("Problem: input should be an integer number!")
        top_n=""
print("\nIMPORTANT is the choice of the receptor. The default option is:",receptor)
receptor_cont=input("Do you want to continue (y/n):")
while receptor_cont=="n":
    receptor=input("What receptor do you want instead?:")
    if os.path.isdir(Path("/home/selina/docking/results")/receptor):
        receptor_cont="y"
    else:
        print("Receptor does not exist in directory: /home/selina/docking/results")
        print("Better options are", [i for i in os.listdir("/home/selina/docking/results") if not i.startswith(".")])
print("\nProblem of using two databanks for ligands is the reoccurrence of the same ligand with different affinity results. The default option for dealing with those is to only keep the result that show the lowest score.")
dublicates=False # default
doubles_cont=input("Do you want to continue (y/n):")
while doubles_cont=="n":
    print("The only other option so far is to simply keep all results included.")
    doubles_cont=input("Do you want to continue (y/n):")
    if doubles_cont=="n":
        print("Too sad. I don't have any other options :(")
    else:
        dublicates=True

################################################ 
# functions
################################################ 
def read_txt_files(file,mode="r"):
    data=open(file, mode=mode).readlines()
    data=[line.replace("\n","") for line in data]
    return data
def write_txt_file(list,file,mode="w"):
    list=[str(line)+"\n" for line in list]
    open(file, mode='w').writelines(list)

def select_row(df,row,values):
    if type(values)==list:
        mask = df[row].isin(values)
    elif type(values) in [int,float,str]:
        mask = df[row] == values
    return df[mask]
    
# combines vina results to csv file (important for ranking etc)
def get_vina_results(receptor,
                     result_version,
                     database,
                     detail=False):
    print("------------------------------")
    # get the ligand results for specified receptor and specified result_version (usually equivilant to database)
    print("Database of ligands:",database)
    path_results=Path("/home/selina/docking/results")
    path_results=path_results/receptor/result_version/"ligand_results"
    ligands=os.listdir(path_results)
    print("Original number of ligands: ", len(ligands))
    # the original files that serve as input for vina
    path_ligandfile=Path("/home/valerie/docking/structures/ligand")/result_version/"pdbqt"
    # when result_version and database variables are not equivilant 
    if not os.path.exists(path_ligandfile):
        path_ligandfile=Path("/home/valerie/docking/structures/ligand")/database/"pdbqt"
    print("ligand files (vina input) from", path_ligandfile,"\n")

    ################# combine vina results
    # define empty variables 
    min_modes={} # dictionary of modes (conformations) with lowest affinity value
    average_values={} # combines all affinity values from given ligand conformations 
    problematic_ligand_files=[] # lists all files where no results from vina are given (problem with input)
    p=0 # counts problematic files
    
    for ligand in ligands:
        try:
            df=pd.read_csv(path_results/ligand/"logfiles"/"results.csv")
            min_affinity=df["affinity"].min()
            mean_affinity=df["affinity"].mean()
            median_affinity=df["affinity"].median()
            min_modes[ligand]=df.iloc[df.index[df["affinity"]==min_affinity][0],:]
            average_values[ligand]={"mean_affinity":mean_affinity, 
                                    "median_affinity":median_affinity}
        except:
            # print(os.listdir(path_results/ligand/"logfiles"))
            problematic_ligand_files.append(str(path_ligandfile/(ligand+".pdbqt")))
            p+=1

    ################# output 
    # dictionary converted to pandas dataframe 
    print("PROBLEM with", p, "ligands (probably missing output files from vina)")
    df_min_mode=pd.DataFrame(min_modes).T
    df_average=pd.DataFrame(average_values).T
    print("Number of ligands with actual vina output: ", df_min_mode.shape[0],"\n")
    
    # list all generated files
    print("created files:")
    path_analysis=Path("/home/selina/docking/results")/receptor/result_version/"analysis"
    if not os.path.isdir(path_analysis): os.makedirs(path_analysis)
    # output: ligand conformations with lowst scores
    new_filename="min_modes.csv"
    print(path_analysis/new_filename)
    df_min_mode.to_csv(path_analysis/new_filename, index=True, index_label="ligand", header=True)
    # output: average values over all conformations of each ligand
    new_filename="average_affinities.csv"
    print(path_analysis/new_filename)
    df_average.to_csv(path_analysis/new_filename, index=True, index_label="ligand", header=True)
    # output: list of ligands where vina output is missing 
    problem_files="problematic_pdbqt_files.txt"
    path=Path("/home/selina/docking/results")/receptor/result_version
    print(path/problem_files)
    write_txt_file(problematic_ligand_files, path/problem_files)

    
    return df_min_mode # type: pd.DataFrame

# convert ZINC and ECBD IDs 
def zinc_and_ecbd_names(ligand_db,ligands,name_file="/home/selina/docking/analysis/ECBD.csv"):
    tab_names=pd.read_csv(name_file)
    tab_names=tab_names[["eos","zinc"]]
    missing_names={}
    
    for ligand in ligands:
        #zinc to ecbd names
        if ligand_db in ["ZINC","zinc"]:
            if ligand in list(tab_names["zinc"]):
                id=tab_names[tab_names["zinc"]==ligand].index[0]
                missing_names[ligand]=tab_names.loc[id,"eos"]  
        #ecbd to zinc names
        elif ligand_db in ["ECBD", "ecbd", "eos"]:
            if ligand in list(tab_names["eos"]):
                id=tab_names[tab_names["eos"]==ligand].index[0]
                missing_names[ligand]=tab_names.loc[id,"zinc"]
                
    return missing_names # type: dictionary

# get dataframe with low scorer of each ligand
def get_low_scorers(receptor,ligand_db,version):

    df_vina=get_vina_results(
        receptor=receptor,
        result_version=version,
        database=ligand_db,
        detail=False
    )

    # add IDs of ZINC or ECBD database (if possible)
    missing_names=zinc_and_ecbd_names(ligand_db,df_vina.index)
    if ligand_db in ["ECBD","ecbd","eos"]:
        df_vina.insert(0, "ZINC", missing_names)
        df_vina=df_vina.reset_index(names="ECBD")
    elif ligand_db in ["ZINC","zinc"]:
        df_vina.insert(0, "ECBD", missing_names)
        df_vina=df_vina.reset_index(names="ZINC")  

    # create csv file with low scorer of each ligand
    path_results=Path("/home/selina/docking/results")
    path=path_results/receptor/version/"analysis"
    if not os.path.isdir(path): os.makedirs(path)
    file=path/("min_modes"+".csv")
    df_vina.to_csv(file, index=False,header=True)

    return df_vina # type: pd.DataFrame

def find_best_ligands(df_vina_results,n):
    top_ligands=df_vina_results.sort_values("affinity")[:n]
    return top_ligands[["ZINC","ECBD","affinity"]]

# finds top scorer for one database (e.g. ZINC)
def get_top_scorers(receptor,ligand_db,top_n,
                   version):
    
    df=get_low_scorers(receptor=receptor,ligand_db=ligand_db,
                       version=version)

    # select top n
    top_ligands=find_best_ligands(df,top_n)
    return top_ligands


def remove_dublicated_rows(df):
    
    dublicates_both=df[df.duplicated(["ZINC","ECBD"],keep=False)]
    no_dublicates_both=df[~df.duplicated(["ZINC","ECBD"],keep=False)]
    dublicates_zinc=df[df.duplicated(["ZINC"],keep=False)]
    dublicates_ecbd=df[df.duplicated(["ECBD"],keep=False)]
    
    max_index_zinc=dublicates_both.groupby(["ZINC"]).apply(
        lambda x: list(x.index[x["affinity"]==max(x["affinity"])]))
    max_index_ecbd=dublicates_both.groupby(["ECBD"]).apply(
        lambda x: list(x.index[x["affinity"]==max(x["affinity"])]))
    
    doubles_dict_zinc={k:v for k,v in Counter(dublicates_zinc["ZINC"]).items() 
                  if not str(k) == "nan"}
    doubles_dict_ecbd={k:v for k,v in Counter(dublicates_ecbd["ECBD"]).items() 
                  if not str(k) == "nan"}
    
    double_with_2_unique_ids_zinc = {k for k,v in doubles_dict_zinc.items() if k not in max_index_zinc.index}
    double_with_2_unique_ids_ecbd = {k for k,v in doubles_dict_ecbd.items() if k not in max_index_ecbd.index}
    
    max_index_zinc_2=select_row(no_dublicates_both,"ZINC",list(double_with_2_unique_ids_zinc)).groupby(["ZINC"]).apply(
        lambda x: list(x.index[x["affinity"]==max(x["affinity"])]))
    max_index_ecbd_2=select_row(no_dublicates_both,"ECBD",list(double_with_2_unique_ids_zinc)).groupby(["ECBD"]).apply(
        lambda x: list(x.index[x["affinity"]==max(x["affinity"])]))

    # combine all and find exact indexes 
    list_min_index=list(set([i[0] for i in max_index_zinc.values] 
                            + [i[0] for i in max_index_ecbd.values]
                            + [i[0] for i in max_index_zinc_2.values]
                            + [i[0] for i in max_index_ecbd_2.values]))

    df_without_doubles=df.drop(list_min_index)
    
    print("We have",len(list_min_index),"ligands with multiple results.")
    
    return df_without_doubles.reset_index(drop=True)
    
# finds top scorer for all database (ZINC & ECBD)
def get_all_top_scorers(receptor,top_n,
                        dublicates=False,
                        zinc={"version":"ZINC",
                              "db":"ZINC"},
                        ecbd={"version":"ECBD",
                              "db":"ECBD"}):
    
    # get complete list of ligands and their resoective lowest affinity score 
    df_zinc=get_low_scorers(receptor=receptor,
                            ligand_db=zinc["db"],
                            version=zinc["version"])
    df_ecbd=get_low_scorers(receptor=receptor,
                            ligand_db=ecbd["db"],
                            version=ecbd["version"])

    # define path for output files 
    path_analysis=Path("/home/selina/docking/analysis")/receptor
    if not os.path.isdir(path_analysis): os.mkdir(path_analysis)

    print("------------------------------")
    print("ANALYSIS esults are in:",path_analysis)
    
    # combine low scorers from ZINC and ECBD 
    combine_all=pd.concat([df_zinc,df_ecbd]).sort_values("affinity").reset_index(drop=True)
    print("All ligands (from ZINC and ECBD) combined:", combine_all.shape[0],"\n")
    all_and_no_dublicates=remove_dublicated_rows(combine_all)
    print("All unique ligands (from ZINC and ECBD) combined:", all_and_no_dublicates.shape[0])
    
    # select top n
    if not dublicates:
        print("Continue with ONLY the low scorers!\n")
        top=all_and_no_dublicates[:top_n]
        name_top="top_unique_"+str(top_n)
    elif dublicates:
        print("Continue to keep all multiple results!\n")
        top=combine_all[:top_n]
        name_top="top_"+str(top_n)
    # create list of ligand molecule sizes (for checking)
    top_size=[]
    for name_z,name_e in zip(top["ZINC"],top["ECBD"]):
        try:
            top_size.append(name_z +" "+ str(get_mol_size(name_z)))
        except:
            top_size.append(name_e +" "+ str(get_mol_size(name_e)))



    # list all generated files
    print("created files:")
    
    # output: csv file with top n ligands
    new_filename=name_top+"_ligands.csv"
    print(path_analysis/new_filename)
    top.to_csv(path_analysis/new_filename, index=False, header=True)
    # output: csv files with all ligands (low scorer of each ligand)
    new_filename="All_ligands.csv"
    print(path_analysis/new_filename)
    combine_all.to_csv(path_analysis/new_filename, index=False, header=True)
    # output: csv files with only unique ligands (no dublicates)
    new_filename="All_unique_ligands.csv"
    print(path_analysis/new_filename)
    all_and_no_dublicates.to_csv(path_analysis/new_filename, index=False, header=True)
    # output: molecule size (x,y,z) for top n ligands
    new_filename=name_top+"_moleculesize.txt"
    print(path_analysis/new_filename)
    write_txt_file(top_size, path_analysis/new_filename)
    
    return top # type: pd.DataFrane

def get_mol_size(name):
    # look for ligand name in original directory (either ZINC or ECBD)
    path=Path("/home/valerie/docking/structures/ligand/ECBD/pdbqt")/(name+".pdbqt")
    if not os.path.exists(path):
        path=Path("/home/valerie/docking/structures/ligand/ZINC/pdbqt")/(name+".pdbqt")
    # compute molecule size
    data=open(path, 'r').readlines()
    coord = [[float(line[31:38]), float(line[39:46]), float(line[47:54])] for line in data if line.split()[0] in ('ATOM', 'HETATM')]
    xcoor, ycoor, zcoor = zip(*coord)
    X,Y,Z = [list(xcoor), list(ycoor), list(zcoor)]
    ranges=[[min(X), max(X)],[min(Y), max(Y)],[min(Z), max(Z)]]
    box = [i for i in [ranges[0][1]-ranges[0][0], ranges[1][1]-ranges[1][0], ranges[2][1]-ranges[2][0]]]
    
    return box # type: list

def get_size_filter(receptor,
                    result_version,
                    file_list="size_filter.txt"):
    
    path=Path("/home/selina/docking/results")/receptor/result_version
    filter_txt=read_txt_files(path/file_list)

    file="conf.txt"
    config=read_txt_files(path/file)
    
    search_box = [float(line[9:]) for line in config if line.startswith("size")]
    search_box_sorted=sorted(search_box)
    
    box_sizes={}
    box_sizes["search box"]={"small":search_box_sorted[0],
                             "medium":search_box_sorted[1],
                             "big":search_box_sorted[2],
                             # "x":search_box[0], "y":search_box[1], "z":search_box[2]
                            }
    for file in filter_txt:
        ligand=re.findall(r"pdbqt/(.*?).pdbqt$",file)[0]
        box=get_mol_size(ligand)
        
        box_sorted=sorted(box)
        box_sizes[ligand]={"small":box_sorted[0],
                           "medium":box_sorted[1],
                           "big":box_sorted[2],
                          }
        
    box_sizes=pd.DataFrame(box_sizes).T
    new_filename=file_list[:-4]+".csv"
    # print(path/new_filename)
    box_sizes.to_csv(path/new_filename, index=True, index_label="ligand", header=True)   
    return box_sizes


    
########################################################

########################################################

# wanted results from analysis: list top scorers
df_top=get_all_top_scorers(receptor,top_n,dublicates=dublicates)

if not dublicates:
    name_top="top_unique_"+str(top_n)
elif dublicates:
    name_top="top_"+str(top_n)
# NEXT: copy of created files into new directory
print("------------------------------")
zinc_names=list(set(df_top["ZINC"]))
ecbd_names=list(set(df_top["ECBD"]))
top_ligand_names=list(set(zinc_names+ecbd_names))

path_results=Path("/home/selina/docking/results")
path_analysis=Path("/home/selina/docking/analysis")
path_valerie_structures=Path("/home/valerie/docking/structures/ligand")

new_dir=Path("/home/selina/docking/top_scorer/top_"+str(top_n)+"_"+receptor)
if not os.path.isdir(new_dir/"all_sdf_files"): os.makedirs(new_dir/"all_sdf_files")
print("COPY of all relevant files to:",new_dir)
print("\t- vina OUTPUT files for each ligand \n",
      "\t- vina conf file \n",
      "\t- vina size filtered ligands \n",
      "\t- ANALYSIS files \n",
      "\t- SDF files (originals)\n",
     )

# go through all ligand names and copy files and directories into new_dir 
# ATTENTION: The databases overlap, 
#            which leads to some ligands having more than one original file and several affinity scores 
i=1
for l in top_ligand_names:
    
    if l in os.listdir(path_results/receptor/"ECBD"/"ligand_results"):
        i+=1
        shutil.copytree(path_results/receptor/"ECBD"/"ligand_results"/l, 
                        new_dir/l,
                        dirs_exist_ok=True)
        shutil.copyfile(path_valerie_structures/"ECBD"/"pilot_library_sdf_files"/(l+".sdf"), 
                        new_dir/"all_sdf_files"/(l+".sdf"))
    elif l in os.listdir(path_results/receptor/"ZINC"/"ligand_results"):
        i+=1
        shutil.copytree(path_results/receptor/"ZINC"/"ligand_results"/l, 
                        new_dir/l,
                        dirs_exist_ok=True)
        shutil.copyfile(path_valerie_structures/"ZINC"/"fda_ligand_sdf_files"/(l+".sdf"), 
                        new_dir/"all_sdf_files"/(l+".sdf"))
        shutil.copyfile(path_valerie_structures/"ZINC"/"fda_with_hydrogens"/(l+".with_hydrogens.sdf"), 
                        new_dir/"all_sdf_files"/(l+"_h.sdf"))

# copy vina output files 
# count of ouput files
print("copy of",i,"output files")
# for ECBD:
file="conf.txt"
shutil.copyfile(path_results/receptor/"ECBD"/file,new_dir/("ECBD_"+file))
file="size_filter.csv"
if not os.path.exists(path_results/receptor/"ECBD"/file): get_size_filter(receptor, "ECBD")
shutil.copyfile(path_results/receptor/"ECBD"/file,new_dir/("ECBD_"+file))
file="problematic_pdbqt_files.txt"
shutil.copyfile(path_results/receptor/"ECBD"/file,new_dir/("ECBD_"+file))
# for ZINC:
file="conf.txt"
shutil.copyfile(path_results/receptor/"ZINC"/file,new_dir/("ZINC_"+file))
file="size_filter.csv"
if not os.path.exists(path_results/receptor/"ZINC"/file): get_size_filter(receptor, "ZINC")
shutil.copyfile(path_results/receptor/"ZINC"/file,new_dir/("ZINC_"+file))
file="problematic_pdbqt_files.txt"
shutil.copyfile(path_results/receptor/"ZINC"/file,new_dir/("ZINC_"+file))

# copy ALL ANALYSIS files 
path_analysis=path_analysis/receptor
new_path_analysis=new_dir/"analysis"
if not os.path.isdir(new_dir/"analysis"): os.makedirs(new_dir/"analysis")
file="All_ligands.csv"
shutil.copyfile(path_analysis/file,new_path_analysis/file)
file="All_unique_ligands.csv"
shutil.copyfile(path_analysis/file,new_path_analysis/file)
file=name_top+"_ligands.csv"
shutil.copyfile(path_analysis/file,new_path_analysis/file)
file=name_top+"_moleculesize.txt"
shutil.copyfile(path_analysis/file,new_path_analysis/file)


