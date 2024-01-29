# Tanimoto
This is how to find analogs of molecules generated with Pocket2Mol. Take your best hits from docking and place them in file SMILES.txt (only their SMILES representations)

## Download database
You can use any database with SMILES format for our project we used ENAMINE Diverse from [ENAMINE](https://enamine.net/compound-collections/real-compounds/real-database-subsets)

## Prepare database
Important part is to convert database from SMILES to fingerprints. For this run

    python generate_fingerprints.py

## Searching database
In order to find analogs run

    python check_similarity

