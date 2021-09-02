import gc
import scanit
import torch
import random
import scanpy as sc
import pandas as pd
import anndata
import numpy as np
from scipy import sparse
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from sklearn.cluster import SpectralClustering, KMeans
import matplotlib.pyplot as plt
import stlearn as st
from pathlib import Path

sp_datadir = './data/slideseq-mouse-cerebellum'

pts = np.loadtxt(sp_datadir+"/positions.csv")
X_sp = sparse.load_npz(sp_datadir+"/expression.npz")
X_sp = X_sp.toarray()
genes_sp = np.loadtxt(sp_datadir+"/genes.txt", dtype = str)
adata = anndata.AnnData(X=X_sp, var=pd.DataFrame(index=genes_sp))
adata.obsm['spatial'] = pts

from somde import SomNode
pts = adata.obsm['spatial']
df_sp = pd.DataFrame(data=adata.X, columns=list(adata.var_names))
som = SomNode(pts, 5)
ndf,ninfo = som.mtx(df_sp.T)
nres = som.norm()
result, SVnum =som.run()
result.to_csv('./data/slideseq-mouse-cerebellum/somde_result.csv')



sp_datadir = './data/slideseq-mouse-hippocampus'

pts = np.loadtxt(sp_datadir+"/positions.csv")
X_sp = sparse.load_npz(sp_datadir+"/expression.npz")
X_sp = X_sp.toarray()
genes_sp = np.loadtxt(sp_datadir+"/genes.txt", dtype = str)
adata = anndata.AnnData(X=X_sp, var=pd.DataFrame(index=genes_sp))
adata.obsm['spatial'] = pts

from somde import SomNode
pts = adata.obsm['spatial']
df_sp = pd.DataFrame(data=adata.X, columns=list(adata.var_names))
som = SomNode(pts, 5)
ndf,ninfo = som.mtx(df_sp.T)
nres = som.norm()
result, SVnum =som.run()
result.to_csv('./data/slideseq-mouse-hippocampus/somde_result.csv')



sp_datadir = './data/slideseq-mouse-olfactory_bulb'

pts = np.loadtxt(sp_datadir+"/positions.csv")
X_sp = sparse.load_npz(sp_datadir+"/expression.npz")
X_sp = X_sp.toarray()
genes_sp = np.loadtxt(sp_datadir+"/genes.txt", dtype = str)
adata = anndata.AnnData(X=X_sp, var=pd.DataFrame(index=genes_sp))
adata.obsm['spatial'] = pts

from somde import SomNode
pts = adata.obsm['spatial']
df_sp = pd.DataFrame(data=adata.X, columns=list(adata.var_names))
som = SomNode(pts, 1)
ndf,ninfo = som.mtx(df_sp.T)
nres = som.norm()
result, SVnum =som.run()
result.to_csv('./data/slideseq-mouse-olfactory_bulb/somde_result.csv')