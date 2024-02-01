import streamlit as st
import pandas as pd
import time
from pathlib import Path
from process_nsp import retrieve_nsp13, nsp13_pdb_to_pdbqt, generate_pocket
from process_ligands import download_ligands, prepare_ligands, convert_ligands_to_pdbqt
from docking import run_docking_for_all_ligands, retrieve_docking_results

# Main function to run the Streamlit app
def main():

    step_times = {}
    st.title("Molecular Docking Application")

    # Sidebar for input parameters
    with st.sidebar:
        st.header("Input Parameters")
        zinc_ids = st.text_input("Enter ZINC IDs (comma-separated)")
        structure_name = st.text_input("Enter structure name")

    # Button to trigger the tasks
    if st.sidebar.button("Run Docking Pipeline"):
        st.write("Running docking pipeline...")
        
        # Step 1: Download molecules
        st.write("Step 1: Downloading molecules...")
        start = time.time()
        retrieve_nsp13(structure_name)
        st.success("Done for nsp13!")
        download_ligands(zinc_ids)
        st.success("Done for ligands!")
        step_times["Download"] = time.time() - start

        # Step 2: Convert molecules
        st.write("Step 2: Converting molecules...")
        start = time.time()
        nsp13_pdb_to_pdbqt(structure_name)
        convert_ligands_to_pdbqt()
        st.success("Done!")
        step_times["Convertion"] = time.time() - start

        # Step 3: Detect pocket
        st.write("Step 3: Detecting pocket...")
        start = time.time()
        generate_pocket()
        step_times["Pocket Detection"] = time.time() - start
        path = "C:/Users/alexi/OneDrive/Bureau/meetuStudent/Meet-U-2023-2024/code/data/pockets/6zsl.cif_predictions.csv"
        pockets = pd.read_csv(path, sep=",")
        st.success("Done!")
        st.dataframe(pockets.head(5))
        clicked = False
        pocket_number = st.selectbox("Which pocket use for docking ?", list(pockets["  rank"].values))
        if pocket_number:
            clicked = True

        
        if clicked:
            # Step 4: Prepare ligands
            st.write("Step 4: Preparing ligands...")
            start = time.time()
            prepare_ligands()
            st.success("Done!")
            step_times["Ligand preparation"] = time.time() - start

            cluster = st.toggle("Run docking on cluster")
            if not cluster:
                # Step 5: Run docking
                start = time.time()
                st.write(f"Step 5: Running docking with pocket number {pocket_number}")
                run_docking_for_all_ligands(pocket_number=pocket_number)
                step_times["Docking"] = time.time() - start

                st.success("Docking pipeline completed successfully!")

                #results = retrieve_docking_results()
                results = {"ZINC001485529916": 5.9, "ZINC001485529917": 6.3}
                results = pd.Series(results)
                results_time = pd.Series(step_times).sort_values()
                st.header("Results:")
                st.bar_chart(results)
                st.header("Time")
                st.bar_chart(results_time)
            
            if cluster:
                cpu = st.number_input("Select number of CPUs")
                nodes = st.number_input("Select number of nodes")
                slurm_file = f"""
##!/bin/bash
#
#SBATCH -N {nodes} # nombre de nœuds
#SBATCH -n {cpu} # nombre de cœurs
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

module load python/3.7
pip install -r requirements.txt

python docking_pipeline.py()
"""
                st.code(slurm_file)

if __name__ == "__main__":
    main()
