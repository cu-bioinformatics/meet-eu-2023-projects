# embedder.py input.fasta output.csv


from bio_embeddings.embed import SeqVecEmbedder
from Bio import SeqIO
import numpy as np
from sys import argv


sequences = [str(record.seq) for record in SeqIO.parse(argv[1], "fasta")]
embedder = SeqVecEmbedder()
embeddings = embedder.embed_many([s for s in sequences])

# reduce_per_protein -> "For a variable size embedding, returns a fixed size embedding encoding all information of a sequence."
reduced_embeddings = np.array([embedder.reduce_per_protein(e) for e in embeddings])
np.savetxt(argv[2], reduced_embeddings, delimiter=",")
