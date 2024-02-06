#!/usr/bin/env python
import argparse
import sys
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

from utils.utils import red_text

SCRIPT_NAME = pathlib.Path(__file__).name

COLORS = [
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
    "#1f78b4",
    "#ffbb78",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=str, help="csv file with ligand embeddings.")
    parser.add_argument("cluster_file", type=str, help="csv file with ligand clusters.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    input_file = pathlib.Path(args.input_file)
    cluster_file = pathlib.Path(args.cluster_file)

    if not input_file.exists():
        print(
            f"{red_text(SCRIPT_NAME)}: {red_text('error')}: {input_file} No such file or directory!"
        )
        return 1
    if input_file.is_dir():
        print(f"{red_text(SCRIPT_NAME)}: error: {input_file} is a directory.")
        return 1
    if not cluster_file.exists():
        print(
            f"{red_text(SCRIPT_NAME)}: {red_text('error')}: {cluster_file} No such file or directory!"
        )
        return 1
    if cluster_file.is_dir():
        print(f"{red_text(SCRIPT_NAME)}: error: {cluster_file} is a directory.")
        return 1

    np.random.seed(11)
    X = pd.read_csv(input_file, index_col=-1)
    clusters = pd.read_csv(cluster_file, index_col=0)

    X = pd.merge(X, clusters, left_index=True, right_index=True)

    _, ax = plt.subplots(figsize=(10, 10))
    pca = PCA(n_components=2)
    x_new = pca.fit_transform(X.drop(columns=["cluster"]).values)

    x_new = pd.DataFrame(x_new, columns=["PCA1", "PCA2"], index=X.index)
    x_new["cluster"] = X["cluster"]
    cl_colors = {i: COLORS[i-1] for i in range(1,13)}
    plt.scatter(
        x_new["PCA1"], x_new["PCA2"], c=x_new["cluster"].map(cl_colors), alpha=0.7
    )
    plt.xlabel("PC1", fontdict={"fontsize": 18})
    plt.ylabel("PC2", fontdict={"fontsize": 18})
    plt.show()
    return 0


if __name__ == "__main__":
    sys.exit(main())
