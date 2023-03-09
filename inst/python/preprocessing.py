"""Preprocessing functions and recipes"""
from typing import Optional, List

import numpy as np
import scanpy as sc
import warnings
from anndata import AnnData
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


class GEXPreprocessor:
    """
    GEX preprocessing pipeline: log1p-transform,
    filter non-expressed cells,
    perform truncated SVD, and row-wise Z-score normalization.
    """

    def __init__(
            self,
            do_log1p: Optional[bool] = True,
            n_components: Optional[int] = 300,
            do_tsvd_before_zscore: Optional[bool] = True,
    ):
        """
        Parameters
        ----------
        do_log1p
            Log1p-transform data?
            Default 'True' expects non-logarithmized data.
        n_components
            Number of components to use for truncated SVD.
        do_tsvd_before_zscore
            Perform truncated SVD before Z-score normalization?
            Default 'True' works better for downstream training and prediction on data from a single dataset.
            Set to 'False' to extract more robust features that work well across datasets.
        """
        self.do_log1p = do_log1p
        self.tsvd = TruncatedSVD(n_components=n_components)
        self.do_tsvd_before_zscore = do_tsvd_before_zscore

    def _transform_before_tsvd(
            self,
            adata: AnnData,
            copy: Optional[bool] = False,
    ) -> AnnData:
        """
        Filter non-expressed cells,
        optionally log1p-transform and z-score normalize.
        Parameters
        ----------
        adata
            Annotated data matrix.
        copy
            Copy adata before transforming?

        Returns
        -------
        Returns the transformed adata.
        """
        if copy:
            adata = adata.copy()
        if self.do_log1p:
            adata.X = sc.pp.log1p(adata.X)

        # filter out non-expressed cells
        with warnings.catch_warnings():
            # ignore UserWarning: Observation names are not unique.
            warnings.filterwarnings("ignore", category=UserWarning)
            sc.pp.filter_cells(adata, min_counts=1)

        if not self.do_tsvd_before_zscore:
            # row-wise Z-score normalization
            adata.X = zscore_normalization(adata.X)
        return adata

    def fit(
            self,
            adata: AnnData,
            copy: Optional[bool] = False,
    ):
        """
        Fit the truncated SVD.
        Default `copy=False` modifies the input AnnData.
        Parameters
        ----------
        adata
            Annotated GEX matrix.
        copy
            Copy adata before fitting?
        """
        self.fit_transform(adata, copy=copy)

    def transform(
            self,
            adata: AnnData,
            copy: Optional[bool] = False,
    ) -> Optional[AnnData]:
        """
        Transform the GEX matrix.
        Parameters
        ----------
        adata
            Annotated GEX matrix.
        copy
            Return a copy if true.
        Returns
        -------
        Returns or updates `adata`, depending on `copy`.
        """
        # filter non-expressed cells, optionally log1p-transform and z-score normalize
        adata = self._transform_before_tsvd(adata, copy=copy)
        X_lowdim = self.tsvd.transform(adata.X)
        if self.do_tsvd_before_zscore:
            X_lowdim = zscore_normalization(X_lowdim)
        adata.obsm["X_pca"] = X_lowdim
        return adata if copy else None

    def fit_transform(
            self,
            adata: AnnData,
            copy: Optional[bool] = False,
    ) -> Optional[AnnData]:
        """
        Fit the truncated SVD and transform the GEX matrix.
        Parameters
        ----------
        adata
            Annotated GEX matrix.
        copy
            Return a copy if true.
        Returns
        -------
        Returns or updates `adata`, depending on `copy`.
        """
        adata = self._transform_before_tsvd(adata, copy=copy)
        X_lowdim = self.tsvd.fit_transform(adata.X)
        if self.do_tsvd_before_zscore:
            X_lowdim = zscore_normalization(X_lowdim)
        adata.obsm["X_pca"] = X_lowdim
        return adata if copy else None
