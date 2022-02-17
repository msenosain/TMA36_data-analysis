import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from pathlib import Path
import requests
import itertools
import re
import sys
from matplotlib.ticker import PercentFormatter
import scanpy as sc
import anndata as ad
import cello
from statannot import add_stat_annotation
from statannotations.Annotator import Annotator
import harmonypy as hm
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_samples, silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import scmer
import pickle


# read adata
results_file = "../../data/TMA36_project/sc_RNAseq/processed/analyzed/all_samples.h5ad"
adata = sc.read(results_file)

# run scmer
model = scmer.UmapL1(
    lasso=3.87e-4,
    ridge=0.0,
    n_pcs=40,
    perplexity=100.0,
    use_beta_in_Q=True,
    n_threads=8,
    pca_seed=2020,
)

# fit model
model.fit(adata.X)

# write model
filename = "../../data/TMA36_project/sc_RNAseq/processed/analyzed/scmer_model"
with open(filename, 'wb') as files:
    pickle.dump(model, files)

