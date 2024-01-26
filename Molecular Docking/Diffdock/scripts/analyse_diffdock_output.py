from pymol import cmd
import os
import numpy as np

dir="/home/valerie/Documents/MeetEU/Diffdock/diffdock_output/"

def get_all_rankfiles_for_ligand(ligand_dir):
    #load and run check_if_ligand_in_binding_site(sel_extent) for all ranks of one ligand
    #test if ligand_dir empty
    sdffiles = []
    #get all sdf_files in output_dir of specific ligand
    for file in os.listdir(ligand_dir):
        if "pdb" not in file:
            sdffiles.append(file)
    #Excess topscorer file needs to be removed
    sdffiles.remove("rank1.sdf")
    #all ranks are safed twice and with two different confidence scores. Only keep file with higher score.
    all_rank_files = []
    for ranknumber in range(1,41):
        temp_list = []
        for filename in sdffiles:
            if "rank"+str(ranknumber)+"_" in filename:
                temp_list.append(filename)
                temp_list.sort()
        all_rank_files.append(temp_list[0])
    return all_rank_files

#check whether all min coordinates (minXYZ) of rank are smaller than min_coordinate of the binding_site_box
# The binding_site box is thereby defined by its center_of_mass. The size of the binding_site box is 30x30x30.
def is_rank_in_box(sel_extent):
    # consensus pocket coordinates for center of mass
    center_x = -14.466000000000001
    center_y = 12.57
    center_z = -73.214
    if sel_extent[0][0]>center_x-15 and sel_extent[0][1]>center_y-15 and sel_extent[0][2]>center_z-15:
        min = True
    else:
        min = False
    #the same as for the min coordinates is repeated for all max coordinates.
    # However, all max coordinates need to be bigger than the max_coordinate of the binding_site_box
    if sel_extent[1][0]<center_x+15 and sel_extent[1][1]<center_y+15 and sel_extent[1][2]<center_z+15:
        max = True
    else:
        max = False

    if min==True and max==True:
        return True
    else:
        return False

#
def check_if_ligand_in_binding_site(ligand_dir):
    all_rank_files = get_all_rankfiles_for_ligand(ligand_dir)
    #load all rank files of one ligand into pymol and get extent values of ligand [[minX,minY,minZ][maxX,maxY,maxZ]]
    #and test if rank is in box
    ligand_in_binding_site = []
    for rank_file in all_rank_files:
        filename_without_ext = os.path.splitext(os.path.basename(rank_file))[0]
        cmd.load(ligand_dir+rank_file)
        sel_extent = cmd.get_extent(str(filename_without_ext))
        in_box = is_rank_in_box(sel_extent)
        ligand_in_binding_site.append(in_box)
    return ligand_in_binding_site

#cound number of complexes for each ligand within binding site
def rank_ligand(dir):
    table_of_diffdock_ranks = []
    for file in os.listdir(dir):
        ligand_id = file
        ligand_dir = dir + str(ligand_id) + "/"
        #check if ligand_dir is empty
        if len(os.listdir(ligand_dir)) == 0:
            print("no files in" + str(ligand_dir))
            continue
        else:
            ranked_ligand = []
            ligand_in_binding_site = check_if_ligand_in_binding_site(ligand_dir)
            count = 0
            for i in ligand_in_binding_site:
                if i == True:
                    count += 1
                if i == False:
                    count += 0
            ranked_ligand.append(ligand_id)
            ranked_ligand.append(count)
            table_of_diffdock_ranks.append(ranked_ligand)
    return table_of_diffdock_ranks

#calculate the percentage of ranks that are within the binding site in relation to percentage of rejected ligands.
#Idealy, top ranked complexes are within binding site.
def calculate_percentage_of_ranks_in_binding_site(ligand_dir):
    ligand_in_binding_site = check_if_ligand_in_binding_site(ligand_dir)
    percentage_of_ranks_in_binding_site = []
    #test if rank1 is in binding site
    for num in range(0,100,10):
        #print(num)
        percentage = num/100
        #print(percentage)
        percentage_of_rejected_ranks = int(40 * percentage)
        #print(percentage_of_rejected_ranks)
        all_ranks_minus_rejected_ranks = ligand_in_binding_site[:40-percentage_of_rejected_ranks]
        #print(len(all_ranks_minus_rejected_ranks))
        count=0
        for i in all_ranks_minus_rejected_ranks:
            if i == True:
                count += 1
            else:
                count += 0
        percentage_of_ranks_in_binding_site.append(round(((count/len(all_ranks_minus_rejected_ranks))*100),2))
    return percentage_of_ranks_in_binding_site

#store results in table
def table_of_percentage(dir):
    table_of_percentage = []
    for file in os.listdir(dir):
        ligand_id = file
        ligand_dir = dir + str(ligand_id) + "/"
        #check if ligand_dir is empty
        if len(os.listdir(ligand_dir)) == 0:
            print("no files in" + str(ligand_dir))
            continue
        else:
            percentage_list = calculate_percentage_of_ranks_in_binding_site(ligand_dir)
            percentage_list.insert(0,ligand_id)
            table_of_percentage.append(percentage_list)
    return table_of_percentage


#table_of_percentage = table_of_percentage(dir)
#table_of_ranked_ligands = rank_ligand(dir)

#export results as csv
"""np.savetxt("table_of_percentage.csv",
           table_of_percentage,
           delimiter =", ",  # Set the delimiter as a comma followed by a space
           fmt ='% s')  # Set the format of the data as string"""

"""np.savetxt("table__ranked_ligands.csv",
           table_of_ranked_ligands,
           delimiter =", ",  # Set the delimiter as a comma followed by a space
           fmt ='% s')  # Set the format of the data as string"""

