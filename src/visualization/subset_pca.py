#!/usr/bin/env python
import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from utils.utils import red_text

SCRIPT_NAME = pathlib.Path(__file__).name


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_file", type=str, help="csv file with sequence embeddings."
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    np.random.seed(11)

    input_file = pathlib.Path(args.input_file)
    if not input_file.exists():
        print(
            f"{red_text(SCRIPT_NAME)}: {red_text('error')}: {input_file} No such file or directory!"
        )
        return 1
    if input_file.is_dir():
        print(f"{red_text(SCRIPT_NAME)}: error: {input_file} is a directory.")
        return 1

    X = pd.read_csv(input_file, header=None, index_col=0)

    _, ax = plt.subplots(figsize=(10, 10))
    pca = PCA(n_components=2)
    x_new = pca.fit_transform(X.values)

    x_new = pd.DataFrame(x_new, columns=["PCA1", "PCA2"], index=X.index)
    print(x_new)
    plt.scatter(x_new["PCA1"], x_new["PCA2"])
    plt.xlabel("PC1", fontdict={"fontsize": 18})
    plt.ylabel("PC2", fontdict={"fontsize": 18})
    plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
