#!/usr/bin/env python
import argparse
import pathlib
import sys

from rdkit import Chem
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
import pandas as pd
import numpy as np


SCRIPT_NAME = pathlib.Path(__file__).name

# Default parameters
FP_SIZE = 2048
RADIUS = 2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file",
        type=str,
        help="input csv file, with ligands. Should contain 'canonical_smiles' and 'molecule_chembl_id' columns.",
    )
    parser.add_argument(
        "-f",
        "--fp-size",
        type=int,
        default=FP_SIZE,
        help="size of the generated fingerprints",
    )
    parser.add_argument(
        "-r",
        "--radius",
        type=int,
        default=RADIUS,
        help="the number of iterations to grow the fingerprint",
    )
    parser.add_argument(
        "-o", "--output", type=str, help="file to which outputs should be written."
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="don't ask for confirmation when overwriting existing files.",
    )
    return parser.parse_args()


def calculate_fingerprints(ligands: pd.DataFrame, fp_generator):
    # radius=2, are roughly equivalent to ECFP4:
    fps = []
    for smile in ligands["canonical_smiles"]:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            fps.append(fp_generator.GetFingerprint(mol))


    # to calculate Tanimoto similarity
    # similarity = DataStructs.TanimotoSimilarity(fps[0], fps[1])

    # as bit vector:
    # print(np.array(fps[0]))
    # indexes of nonzero values:
    # print(np.nonzero(fps[0]))

    # write as bit vectors: 0:2048 columns are bit representation, last column is molecule_chembl_id
    df = pd.DataFrame(np.array(fps))
    df.insert(0, "molecule_chembl_id", ligands["molecule_chembl_id"])
    return df


def main() -> int:
    args = parse_args()
    input_file = pathlib.Path(args.input_file)

    if not input_file.exists():
        print(f"{SCRIPT_NAME}: error: {input_file}: No such file or directory")
        return 1
    if input_file.is_dir():
        print(f"{SCRIPT_NAME}: error: {input_file} is a directory.")
        return 1

    ligands = pd.read_csv(input_file)
    fpgen = GetMorganGenerator(radius=args.radius, fpSize=args.fp_size)

    df = calculate_fingerprints(ligands, fpgen)

    if args.output:
        output_file = pathlib.Path(args.output)

        if output_file.is_dir():
            print(f"{SCRIPT_NAME}: error: {output_file} is a directory.")
            return 1
        if (not args.force) and output_file.exists():
            print(f"[warning] File {output_file} exists.", end=" ")
            choice = input("Do you want to overwrite it? [Y/n]: ")

            if choice.lower() not in ["y", "yes"]:
                print(df.to_csv(index=False))
                return 0

        df.to_csv(output_file, index=False)
        return 0

    print(df.to_csv(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
