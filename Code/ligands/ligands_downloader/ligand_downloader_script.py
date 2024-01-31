import pandas as pd
import os

# Downloading of the CSV file
csv_file = './ligands/ligands_downloader/csv_files/substances.csv'

#Dataframing w pandas
df = pd.read_csv(csv_file,sep = ';')
#print(df)

# Creation of a folder for the ligands
ligands_folder = './ligands/zinc_database'
os.makedirs(ligands_folder, exist_ok=True)

# Downloading of the SDF Files from the ZINC databse
for zinc_id in df['zinc_id']:
    sdf_url = f'https://zinc.docking.org/substances/{zinc_id}.sdf'
    command = f'wget {sdf_url} -O {ligands_folder}/{zinc_id}.sdf'
    result = os.system(command)
    if result != 0:
        print(f"Échec du téléchargement pour ZINC ID {zinc_id}")
    else:
        print(f"Téléchargement réussi pour ZINC ID {zinc_id}")

