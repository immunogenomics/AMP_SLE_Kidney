import anndata as ad
from starcat import starCAT
import pandas as pd
from scipy.io import mmread

X = mmread("onlyT_cells.mtx").tocsr()  # or .tocsc()

obs = pd.read_csv("onlyT_rows.txt", header=None, names=["cell_id"])
var = pd.read_csv("onlyT_cols.txt", header=None, names=["gene_symbol"])

adata = ad.AnnData(
    X=X,
    obs=obs.set_index("cell_id"),
    var=var.set_index("gene_symbol")
)

tcat = starCAT(reference='TCAT.V1')
usage, scores = tcat.fit_transform(adata.T)

usage.to_csv("tcat_usage.csv")
scores.to_csv("tcat_scores.csv")