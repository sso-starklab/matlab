
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import pandas as pd
import xgboost as xgb
from xgboost import plot_tree, DMatrix
import joblib


def predict(feat_mat):
    xg_class = joblib.load('XGens')
    a = xg_class.predict_proba(feat_mat)
    return a



