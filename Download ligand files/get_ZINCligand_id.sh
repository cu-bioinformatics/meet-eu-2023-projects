#directory in which txt_file from database is supposed to be safed in
dir="/home/valerie/docking/structures/ligand/ZINC"

#get txt_file from ZINC database BUT ONLY FDA approved ones
wget -O $dir/ZINC_fda.txt "https://zinc15.docking.org/substances/subsets/fda.txt?count=all" -P $dir

# downloaded file (fda.txt) is looped through, the ZINC-id is extracted
filename=$dir/ZINC_fda.txt

while IFS= read -r line
do
   echo "${line:0:16}"
done < "$filename"

#IMPORTANT: When running the script make sure the output is stored in a file called fda_idlsit-txt as shown below!
#bash get_ZINCligand_id.sh > /home/valerie/docking/structures/ligand/ZINC/fda_idlist.txt
