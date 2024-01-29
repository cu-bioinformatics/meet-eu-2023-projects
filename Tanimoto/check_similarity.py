# Piotr Trzaskowski
# Functionality:
#     Fingerprint Comparison: The script compares chemical fingerprints, represented as bit vectors, of molecules. It uses the Tanimoto coefficient, a common similarity metric in cheminformatics, to measure the similarity between pairs of fingerprints.
#     Efficient Large-Scale Processing: To handle large databases, the script employs multiprocessing, enabling it to process large fingerprint files in parallel. This is critical for working with extensive databases like the Enamine REAL database, which contains millions of chemical structures.
#     Chunk-Based File Reading: It reads large fingerprint files in manageable chunks (default set to 10MB), ensuring that the memory usage is optimized even when working with very large files.
#     Customizable Parameters: Key parameters like the similarity threshold, chunk size, and number of worker processes are customizable, allowing users to tailor the script to their specific hardware and requirements.
#     Progress Monitoring: The script provides feedback on the progress of the processing, which is crucial for long-running tasks.
#     Output Generation: Results, including pairs of similar molecules and their similarity scores, are written to a CSV file for further analysis.
# Use Case:
# The primary use case of this script is in drug discovery and chemical research, where identifying similar molecules from a large database is a common task. For example, a researcher might use this script to find molecules in the Enamine REAL database that are similar to a set of known active compounds. This could aid in the identification of new drug candidates or in the study of structure-activity relationships.

import multiprocessing
from rdkit.DataStructs.cDataStructs import ExplicitBitVect
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
from datetime import datetime

# the input files should contain SMILES strings and base64 strings of fingerprints separated by a comma and a newline character
# the encoding should be utf-8
your_fingerprint_path = 'SMILES_fingerprints.csv'
database_fingerprint_path = 'Enamine_Diverse_REAL_drug-like_48.2M_cxsmiles_fingerprints.csv'

least_acceptable_similarity = 0.45 # Two structures are usually considered similar if T > 0.85
number_of_bits_per_fingerprint = 2048 # must be the same as in generate_fingerprints.py

chunk_size = 10*1024*1024  # 10MB 
max_workers = 8 # 8 cores powerfull processor

out_file_name = f'similarities_{datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}.csv'

def chunk_to_fingerprints(chunk):
    pairs = []

    # Process each line in the chunk
    for line in chunk.split(b'\n'):
        if line:
            smiles_bytes, fp_bytes = line.strip().split(b',', 1)
            smiles = smiles_bytes.decode('utf-8')
            base64_str = fp_bytes.decode('utf-8')

            # Convert base64 string to ExplicitBitVect
            fp_vect = ExplicitBitVect(number_of_bits_per_fingerprint)
            fp_vect.FromBase64(base64_str)

            pairs.append((smiles, fp_vect))

    return pairs


def read_chunk(file_path): 
    with open(file_path, 'rb') as f:
        chunk = b''  # Initialize an empty bytes object for the chunk
        while True:
            # Read the next chunk
            new_data = f.read(chunk_size)
            if not new_data:
                break
            chunk += new_data

            # Ensure the chunk ends at a line break
            while not chunk.endswith(b'\n'):
                additional_data = f.read(1)
                if not additional_data:
                    break
                chunk += additional_data

            yield chunk
            # Reset chunk for the next iteration
            chunk = b''


def worker(chunk, out_file_name, lock ,lock_print, processed_lines, total_line_count, semaphore):
    base_pairs = chunk_to_fingerprints(chunk)

    buffer = bytearray()

    for chunk in read_chunk(your_fingerprint_path):
        your_pairs = chunk_to_fingerprints(chunk)
        for base_smiles, base_fp in base_pairs:
            for your_smiles, your_ftp in your_pairs:
                tanimoto_coeff = TanimotoSimilarity(your_ftp, base_fp)
                if tanimoto_coeff >= least_acceptable_similarity:
                    buffer.extend(base_smiles.encode('utf-8'))
                    buffer.extend(b',')
                    buffer.extend(your_smiles.encode('utf-8'))
                    buffer.extend(b',')
                    buffer.extend(str(tanimoto_coeff).encode('utf-8'))
                    buffer.extend(b'\n')


    if buffer:
        with lock:
            with open(out_file_name, 'ab') as f:
                f.write(buffer)
        buffer.clear()

    with lock_print:
        processed_lines.value += len(base_pairs)
        current_percentage = (processed_lines.value / total_line_count) * 100
        print(f"Processed: {current_percentage:.2f}%")

    semaphore.release()


if __name__ == "__main__":
    print("calculating lines")

    # Get number of lines in the input file
    total_line_count = 0
    with open(database_fingerprint_path, 'rb') as f:
        while chunk := f.read(chunk_size):  # read in chunks of 1MB
            total_line_count += chunk.count(b'\n')

    print("lines: ", total_line_count)
    
    # for multiprocessing
    processed_lines = multiprocessing.Value('q', 0)
    lock = multiprocessing.Lock()
    lock_print = multiprocessing.Lock()

    semaphore = multiprocessing.Semaphore(max_workers)

    print("processing started")

    p = None

    for chunk in read_chunk(database_fingerprint_path):
        semaphore.acquire()
        p = multiprocessing.Process(target=worker, args=(chunk, out_file_name, lock, lock_print, processed_lines, total_line_count, semaphore))
        p.daemon = True
        p.start()

    for i in range(max_workers):
        semaphore.acquire()

    print("Processing complete.")
