import numpy as np
import pandas as pd
import xgboost as xgb
from xgboost import plot_tree, DMatrix
import joblib
import time
import math

from prepare_features import prepare_features
from uptade_all import update_all


def predict(feat_mat):
    xg_class = joblib.load('XGens')
    a = xg_class.predict_proba(feat_mat)
    return a


def load_data(filebase, shank_num):
    f = open(filebase+'/info')
    info = f.read()
    lst_info = info.split('\n')
    session_name = lst_info[0]

    clu = np.load(filebase + '/' + session_name + '.clu.' + shank_num + '.npy')
    res = np.load(filebase + '/' + session_name + '.res.' + shank_num + '.npy')
    cc = np.load(filebase + '/' + session_name + '.cc.' + shank_num + '.npy')
    mspk = np.load(filebase + '/' + session_name + '.mspk.' + shank_num + '.npy')
    sspk = np.load(filebase + '/' + session_name + '.sspk.' + shank_num + '.npy')
    nspk_vec = np.load(filebase + '/' + session_name + '.nspk_vec.' + shank_num + '.npy')
    nspk_vec = np.squeeze(nspk_vec)
    return clu, res, cc, mspk, sspk, nspk_vec


def main_loop(filebase, shank_num):
    clu, res, cc, mean_spk, std_spk, nspk_vec = load_data(filebase, shank_num)
    n_original = len(clu)

    idx_ignored = np.where((clu == 0) | (clu == 1))
    clu = np.delete(clu, idx_ignored[0], 0)

    thresholds = [0.98, 0.9, 0.8, 0.7, 0.6, 0.5]
    unique_clu = np.unique(clu)
    n = len(unique_clu)
    c = 0
    max_pred = 1
    pairs = []
    while max_pred > 0.5:
        unique_clu = np.unique(clu)
        n = len(unique_clu)
        pairs = []
        feat_mat = np.empty((1, 422))
        for j in range(n):
            for m in range(n):
                if j != m:
                    x = prepare_features(j, m, clu, mean_spk, std_spk, cc, unique_clu)
                    feat_mat = np.append(feat_mat, np.array([x.T]), axis=0)
                    pairs.append((j, m))

        feat_mat = np.delete(feat_mat, 0, axis=0)
        a = predict(feat_mat)
        a = np.round(a, 3)
        max_pred = max(a[:, 1])
        idx = np.where(a[:, 1] == max_pred)
        idx = idx[0][0]
       # u1 = math.floor(idx/n)
       # u2 = (idx % n)
        u1 = pairs[idx][0]
        u2 = pairs[idx][1]
        u1 = unique_clu[u1]
        u2 = unique_clu[u2]

        group_units = np.array([u1, u2])
        clu, mean_spk, std_spk, nspk_vec, cc = update_all(cc, group_units, clu, mean_spk, std_spk, nspk_vec, [0, 1])


    new_clu = clu_organize(filebase, clu, idx_ignored, n_original, shank_num)
    return


def clu_organize(filebase, sorted_clu, idx_ignored, n_original,shank_num):
    new_clu = np.ones((n_original, 1))
    new_clu[idx_ignored] = 0
    c = 0
    for i in range(n_original):
        if new_clu[i][0] == 1:
            new_clu[i][0] = sorted_clu[c]
            c += 1
    save_clu(filebase, new_clu, shank_num)
    return new_clu


def save_clu(filebase, new_clu, shank_num):
    f = open(filebase + '/info')
    info = f.read()
    lst_info = info.split('\n')
    session_name = lst_info[0]

    fname = filebase + "/" + session_name + ".clu." + shank_num + "_2"
    np.save(fname, new_clu)
    return


filebase = '/home/tali/matlab/AUSS_python/npy_files'
f = open(filebase + '/info')
info = f.read()
lst_info = info.split('\n')
n_shanks = int(lst_info[1])
for i in range(n_shanks):
    shank_num = i+1
    shank_num = str(shank_num)
    main_loop(filebase, shank_num)

