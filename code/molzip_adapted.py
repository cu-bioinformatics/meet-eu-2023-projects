import numpy as np
import gzip
from sklearn.decomposition import PCA

"""
This code has been adapted from https://github.com/daenuprobst/molzip.git
"""

def regression(val_smiles, train_smiles, train_val, k=1):
    "Perform the regression"
    # compress first the training data to save some computational time
    ctrain_smiles = [len(gzip.compress(smile.encode())) for smile in train_smiles]
    pred_l = []
    for smile in val_smiles:
        pred_l.append(regression_one(smile, train_smiles, ctrain_smiles, train_val, k))
    return np.array(pred_l)


def projection(val_smiles, train_smiles=None, n_comp=2):
    "Perform a PCA on using the distance matrix"
    if train_smiles is None:
        train_smiles = val_smiles
    # compress first the training data to save some computational time
    ctrain_smiles = [len(gzip.compress(smile.encode())) for smile in train_smiles]
    dist_mat = []
    for smile in val_smiles:
        dist_mat.append(regression_one(smile, train_smiles, ctrain_smiles, [], k=10, proj=True))
    dist_mat = np.array(dist_mat)
    pca = PCA(n_components=n_comp)
    return pca.fit_transform(dist_mat)


def regression_one(x1, X_train, cX_train, y_train, k, proj=False):
    """
    Compute the distance of x1 with the training data

    x1 = target smiles
    X_train = list of training smiles
    cX_train = list of compressed training smiles
    y_train = list of values to fit
    """
    y_train = np.array(y_train)
    # size of the compressed input smile x1
    Cx1 = len(gzip.compress(x1.encode()))
    distance_from_x1 = []

    # compute distances from the training data
    for x2, Cx2 in zip(X_train, cX_train):
        Cx1x2 = len(gzip.compress((x1 + x2).encode()))
        ncd = (Cx1x2 - min(Cx1, Cx2)) / max(Cx1, Cx2)
        distance_from_x1.append(ncd)

    if proj:
        return distance_from_x1
    else:
        distance_from_x1 = np.array(distance_from_x1)

        # find the training smiles that are the closest to the input x1
        sorted_idx = np.argsort(distance_from_x1)
        top_k_values = y_train[sorted_idx[:k]]  # take their values
        top_k_dists = distance_from_x1[sorted_idx[:k]]  # save their disatnce

        # now compute a distance average
        sum_dist = top_k_dists.sum()

        return np.matmul(top_k_dists, top_k_values) * (1./sum_dist)

