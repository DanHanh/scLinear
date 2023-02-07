"""Preprocessing functions and recipes"""
import scanpy as sc
import numpy as np
from anndata import AnnData
from typing import Optional, List
from sklearn.decomposition import TruncatedSVD


def filter_to_common_genes(
        adatas: List[AnnData],
        copy: Optional[bool] = False,
) -> List[AnnData]:
    """
    Filter all AnnDatas to the intersection of their genes.
    Parameters
    ----------
    adatas
        List of AnnDatas.
    copy
        Return a copy if true.
    """
    if copy:
        adatas = [adata.copy() for adata in adatas]
    common_genes = list(set.intersection(*[set(adata.var_names) for adata in adatas]))
    for adata in adatas:
        adata._inplace_subset_var(np.isin(adata.var_names, common_genes))
    return adatas


def recipe_lance22(
        adata: AnnData,
        n_tsvd_comp: Optional[int] = 300,
        log: Optional[bool] = True,
        copy: Optional[bool] = False,
) -> Optional[AnnData]:
    """
    Normalization and filtering as of [Lance22]_.
    Default `log=True` expects non-logarithmized data.
    If using logarithmized data, pass `log=False`.
    Parameters
    ----------
    adata
        Annotated data matrix.
    n_tsvd_comp
        Number of components to use for truncated SVD.
    log
        Logarithmize data?
    copy
        Return a copy if true.
    """
    if copy:
        adata = adata.copy()
    # filter out duplicate variables
    adata.var_names_make_unique()
    # filter out non-expressed genes and cells
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_cells(adata, min_counts=1)
    if log:
        # log(1+x) transform
        adata.X = sc.pp.log1p(adata.X)
    # truncated SVD
    tsvd = TruncatedSVD(n_components=300)
    X_lowdim = tsvd.fit_transform(adata.X.toarray())
    # scale data to unit variance and shift to zero mean (z-score normalization)
    X_lowdim_norm = zscore_normalization(X_lowdim)
    adata.obsm["X_pca"] = X_lowdim_norm
    return adata if copy else None


def zscore_normalization(
        X: np.ndarray
) -> np.ndarray:
    """
    Row-wise Z-score normalization.
    Parameters
    ----------
    X
        Data matrix.
    """
    X_sd = np.std(X, axis=1).reshape(-1, 1)
    X_sd[X_sd == 0] = 1
    X_normalized = (X - np.mean(X, axis=1).reshape(-1, 1)) / X_sd
    return X_normalized.astype(np.float32)
