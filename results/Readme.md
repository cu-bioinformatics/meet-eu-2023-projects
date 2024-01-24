# Docking site classification results

## Guide to result file names

- `filtr_out_[peptide_alias]_[tool_alias].txt` - text files with results of the docked site classification
- `binding_pockets_1.txt` - input file for binding site classification, contains residues involved in the formation of the binding pockets

### Peptide aliases

- `P` - peptides from the literature review
- `K` - peptides from structural search of the PDB database
- `G` - peptides generated using protein language model

If results contain more than one set of peptides the aliases will be concatenated.

If the results contain a subset of peptides, alias will be followed by the range of indices of those peptides (e.g. P1-18 for peptides from P1 to P18).

### Tool aliases

- `a` - AlphaFold2
- `c` - CABS-dock
- `h` - HPepDock
