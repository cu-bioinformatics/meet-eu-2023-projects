# get sdf file for pilot libary from ECBD
dir="/home/valerie/docking/structures/ligand/ECBD"


wget --no-check-certificate "https://ecbd.eu/static/core/compounds/pilot_library.sdf" -P $dir

# -m: split bigsdf file into single ones, -h: add hydrogens (i.e. make them all implicit), --gen3d creates 3d structure
obabel -isdf $dir/pilot_library.sdf -osdf -O $dir/ECBD_sdf_files/.sdf --split -h --gen3d

