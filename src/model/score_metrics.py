"""The script is copied from https://github.com/HUBioDataLab/ProtBENCH/blob/main/scripts/score_metrics.py
"""


"""The functions used for the calculation of score metrics below except from 
pchembl_multiclass and multiclass_mcc are retrieved from the source code of the
IDG-DREAM Drug-Kinase Binding prediction Challenge, which is now published in 
Nature Communications (https://doi.org/10.1038/s41467-021-23165-1)"""


# -*- coding: utf-8 -*-
import numpy as np
import copy

from math import sqrt
from scipy import stats

from sklearn.metrics import f1_score, matthews_corrcoef, accuracy_score
from sklearn.preprocessing import binarize


def rmse(y, y_pred):
    rmse = sqrt(((y - y_pred) ** 2).mean(axis=0))
    return rmse


def pearson(y, y_pred):
    pr = np.corrcoef(y, y_pred)[0, 1]
    return pr


def spearman(y, y_pred):
    sp = stats.spearmanr(y, y_pred)[0]
    return sp


def accuracy(y, y_pred, pchembl):
    y_binary = copy.deepcopy(y)
    y_binary = binarize(y_binary.reshape(1, -1), threshold=pchembl, copy=False)[0]
    y_pred_binary = copy.deepcopy(y_pred)
    y_pred_binary = binarize(
        y_pred_binary.reshape(1, -1), threshold=pchembl, copy=False
    )[0]

    acc = accuracy_score(y_binary, y_pred_binary)
    return acc


def f1(y, y_pred, pchembl):
    y_binary = copy.deepcopy(y)
    y_binary = binarize(y_binary.reshape(1, -1), threshold=pchembl, copy=False)[0]
    y_pred_binary = copy.deepcopy(y_pred)
    y_pred_binary = binarize(
        y_pred_binary.reshape(1, -1), threshold=pchembl, copy=False
    )[0]

    f1 = f1_score(y_binary, y_pred_binary)
    return f1


def mcc(y, y_pred, pchembl):
    y_binary = copy.deepcopy(y)
    y_binary = binarize(y_binary.reshape(1, -1), threshold=pchembl, copy=False)[0]
    y_pred_binary = copy.deepcopy(y_pred)
    y_pred_binary = binarize(
        y_pred_binary.reshape(1, -1), threshold=pchembl, copy=False
    )[0]

    mcc = matthews_corrcoef(y_binary, y_pred_binary)
    return mcc


def pchembl_multiclass(pchembl_list):
    pc_class_list = []
    for pc in pchembl_list:
        if pc < 5.0:
            pc_class = "<5.0"
        elif pc >= 5.0 and pc < 5.5:
            pc_class = "5.0-5.5"
        elif pc >= 5.5 and pc < 6.0:
            pc_class = "5.5-6.0"
        elif pc >= 6.0 and pc < 6.5:
            pc_class = "6.0-6.5"
        elif pc >= 6.5 and pc < 7.0:
            pc_class = "6.5-7.0"
        else: # pc >= 7.0:
            pc_class = "7.0>="
        pc_class_list.append(pc_class)
    return pc_class_list


def multiclass_mcc(y, y_pred):
    mc_y = pchembl_multiclass(y)
    mc_y_pred = pchembl_multiclass(y_pred)

    mc_mcc = matthews_corrcoef(mc_y, mc_y_pred)
    return mc_mcc
