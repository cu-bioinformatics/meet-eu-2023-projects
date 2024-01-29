# This Python script is designed for processing large datasets of chemical compounds represented in SMILES (Simplified Molecular Input Line Entry System) notation. Specifically, it's tailored for generating molecular fingerprints from SMILES strings. Molecular fingerprints are a way of encoding the structure of a molecule in a format that's easy for computer algorithms to handle, often used in cheminformatics for tasks like similarity searching or machine learning.

# Key Features of the Script:

#     Efficient Handling of Large Files: The script is optimized to process large files efficiently, by reading the file in chunks, thereby reducing memory usage.
#     Parallel Processing: Utilizes multiprocessing for parallel processing, significantly speeding up the computation by generating fingerprints concurrently.
#     Molecular Fingerprint Generation: Converts SMILES to molecular fingerprints using RDKit, a widely used toolkit in cheminformatics.
#     Progress Tracking: Includes a mechanism to track and print the progress of processing.

# Use Case:
# This script is particularly useful in drug discovery and chemical informatics, where handling large databases of chemical compounds is common. Researchers can use it to preprocess compound databases, converting them into a format suitable for various computational chemistry tasks like virtual screening, similarity searches, or input for machine learning models.

# Author: Piotr Trzaskowski

# Note: This script assumes familiarity with RDKit and multiprocessing in Python. Users should have a basic understanding of cheminformatics to effectively utilize this tool.

from rdkit import Chem
from rdkit.Chem import AllChem
import multiprocessing

# the input file should contain SMILES strings separated by a newline character
input_filename = 'Enamine_Diverse_REAL_drug-like_48.2M_cxsmiles.cxsmiles'
input_filename = 'SMILES.txt'
# Removing the extension and adding "_fingerprints.csv"
output_filename = '.'.join(input_filename.split('.')[:-1]) + "_fingerprints.csv"

chunk_size=10*1024*1024 # 10MB 
max_workers = 8 # 8 cores powerfull processor

number_of_bits_per_fingerprint = 2048

def generate_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, number_of_bits_per_fingerprint)

def read_chunk(input_file):
    chunk = []
    current_size = 0
    numer_of_lines = 0
    with open(input_file, 'r') as file:
        for line in file:
            numer_of_lines += 1
            smile = line.strip().split()[0]
            chunk.append(smile)
            current_size += len(smile)
            if current_size >= chunk_size:
                yield chunk, numer_of_lines
                
                chunk = []
                current_size = 0
                numer_of_lines = 0
        if chunk:
            yield chunk, numer_of_lines


def process_chunk(chunk, output_file, lock, lock_print, numer_of_lines, processed_lines, total_line_count, semaphore):
    buffer = bytearray() # it's just a char[] wrapper from C language

    for smiles in chunk:
        fp = generate_fingerprint(smiles)
        if fp is not None:
            buffer.extend(smiles.encode('utf-8'))
            buffer.extend(b',')
            buffer.extend(fp.ToBase64().encode('utf-8'))
            buffer.extend(b'\n')

    with lock:
        with open(output_file, 'ab') as f:
            f.write(buffer)
    
    buffer.clear()

    with lock_print:
        processed_lines.value += numer_of_lines
        current_percentage = (processed_lines.value / total_line_count) * 100
        print(f"Processed: {current_percentage:.2f}%")

    semaphore.release()


if __name__ == "__main__":

    # Clear or create the output file
    open(output_filename, 'w').close()

    print("calculating lines")
    # Get number of lines in the input file
    total_line_count = 1
    with open(input_filename, 'r') as f:
        total_line_count = max(sum(1 for line in f), 1)

    print("lines: ", total_line_count)

    # for multiprocessing
    processed_lines = multiprocessing.Value('q', 0)
    lock = multiprocessing.Lock()
    lock_print = multiprocessing.Lock()

    semaphore = multiprocessing.Semaphore(max_workers)  # Semaphore to control chunk reading

    p = None

    print("processing started")

    for chunk, numer_of_lines in read_chunk(input_filename):
        semaphore.acquire()
        p = multiprocessing.Process(target=process_chunk, args=(chunk, output_filename, lock, lock_print, numer_of_lines, processed_lines, total_line_count, semaphore))
        p.daemon = True
        p.start()

    for i in range(max_workers):
        semaphore.acquire()
        
    print("Processing complete.")
