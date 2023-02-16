"""High level functions for running the model prediction pipeline."""
from typing import Optional

import anndata as ad
import numpy as np
from sklearn.linear_model import LinearRegression


class ADTPredictor:
    """Wrapper for running the GEX to ADT prediction pipeline."""

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
            Logarithmize data?
            Default 'True' expects non-logarithmized data.
        n_components
            Number of components to use for truncated SVD.
        do_tsvd_before_zscore
            Perform truncated SVD before Z-score normalization?
            Default 'True' works better for downstream training and prediction on data from a single dataset.
            Set to 'False' to extract more robust features that work well across datasets.

        """
        self.gex_preprocessor = GEXPreprocessor(
            do_log1p=do_log1p,
            n_components=n_components,
            do_tsvd_before_zscore=do_tsvd_before_zscore,
        )
        self.model = LinearRegression()

    def fit(
            self,
            gex_train: np.ndarray,
            adt_train: np.ndarray,
            gex_test: Optional[np.ndarray] = None,
    ):
        """
        Fit the GEX preprocessing and the GEX to ADT model to the training data.
        gex_test is optional and is used for transductive preprocessing,
        i.e. the truncated SVD is fit on the union of the training and test data.
        Parameters
        ----------
        gex_train
            Training GEX data.
        adt_train
            Training ADT data.
        gex_test
            Optional test GEX data for transductive preprocessing.
        """
        # Preprocess GEX data
        gex_train = ad.AnnData(gex_train)
        if gex_test is not None:
            gex_test = ad.AnnData(gex_test)
            gex_train_test = ad.concat((gex_train, gex_test), join='outer')
            self.gex_preprocessor.fit_transform(gex_train_test)
            X_train = gex_train_test.obsm['X_pca'][:gex_train.shape[0]]
        else:
            self.gex_preprocessor.fit_transform(gex_train)
            X_train = gex_train.obsm['X_pca']

        # Fit the model
        self.model.fit(X_train, adt_train)

    def predict(
            self,
            gex_test: np.ndarray,
    ) -> np.ndarray:
        """
        Predict ADT from GEX.
        Parameters
        ----------
        gex_test
            Test GEX matrix.
        Returns
        -------
        adt_pred
            Predicted ADT matrix.
        """
        # Preprocess GEX data
        gex_test = ad.AnnData(gex_test)
        self.gex_preprocessor.transform(gex_test)
        X_test = gex_test.obsm['X_pca']

        # Predict ADT data
        adt_pred = self.model.predict(X_test)
        # Clip negative values
        adt_pred = np.clip(adt_pred, a_min=0, a_max=None)

        return adt_pred
