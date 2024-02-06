
# directory to file where fda_idlist.txt is situated
dir="/home/valerie/docking/structures/ligand/ZINC"

#downloading one sdf file containing all substances did not work.
#sdfurl="https://zinc15.docking.org/substances/subsets/fda.sdf?count=all"

#instead downloading every structur individually:
url="https://zinc15.docking.org/substances/"

# id file contains all ligandids we want to download
id_file="fda_idlist.txt"

while IFS="" read -r ligandid
do
	wget $url$ligandid.sdf -P $dir/zinc_sdf_files/
	echo "downloaded $ligandid from ZINC into $dir/zinc_sdf_files"
done < $dir/$id_file



