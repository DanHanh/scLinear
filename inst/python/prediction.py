"""High level functions for running the model prediction pipeline."""
from typing import Optional, Tuple

import anndata as ad
import numpy as np
from numpy import ndarray
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import LinearRegression
import warnings

from preprocessing import GEXPreprocessor

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
        self.gex_names = None
        self.adt_names = None

    def fit(
            self,
            gex_train: np.ndarray,
            adt_train: np.ndarray,
            gex_test: Optional[np.ndarray] = None,
            gex_names: Optional[np.ndarray] = None,
            adt_names: Optional[np.ndarray] = None,
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
        gex_names
            Optional GEX gene names.
        adt_names
            Optional ADT protein names.
        """
        # If GEX or ADT names are provided, save them
        # ADT names will be returned in the prediction
        if gex_names is not None:
            self.gex_names = gex_names
        if adt_names is not None:
            self.adt_names = adt_names
        # Preprocess GEX data
        gex_train = ad.AnnData(gex_train, dtype=gex_train.dtype)

        if gex_test is not None:
            gex_test = ad.AnnData(gex_test, dtype=gex_test.dtype)
            with warnings.catch_warnings():
                # ignore UserWarning: Observation names are not unique.
                warnings.filterwarnings("ignore", category=UserWarning)
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
            gex_names: Optional[np.ndarray] = None,
    ) -> Tuple[ndarray, ndarray]:
        """
        Predict ADT from GEX.
        Parameters
        ----------
        gex_test
            Test GEX matrix.
        gex_names
            Optional GEX gene names of the columns of `gex_test`.
            If provided, the function will check if `gex_names` matches `self.gex_names` used for training,
            and will use only the columns that match.
            The names in `gex_names` that do are not in `self.gex_names` will be ignored,
            and the columns of `gex_test` that do not have a matching name in `gex_names` will be set to 0.
        Returns
        -------
        adt_pred
            Predicted ADT matrix.
        adt_names
            ADT protein names of the prediction.
        """
        # If GEX names are provided, check if they match the training names
        if gex_names is not None:
            if self.gex_names is None:
                raise ValueError(
                    'GEX names were not provided during training. '
                    'Please provide GEX names during training or prediction '
                    'if you want to match gene names by providing gex_names as an argument.'
                )
            if not np.array_equal(gex_names, self.gex_names):
                # Discard the genes that are not in the training data
                # Set the genes that are in the test data but not in the training data to 0
                # And sort the columns to match the training data
                selfgex2idx = dict()
                for i, g in enumerate(self.gex_names):
                    selfgex2idx[g] = i
                gex_test_new = np.zeros((gex_test.shape[0], len(self.gex_names)))
                for i, g in enumerate(gex_names):
                    if g in selfgex2idx:
                        gex_test_new[:, selfgex2idx[g]] = gex_test[:, i]
                gex_test = gex_test_new
        # Preprocess GEX data
        gex_test = ad.AnnData(gex_test, dtype=gex_test.dtype)
        self.gex_preprocessor.transform(gex_test)
        X_test = gex_test.obsm['X_pca']

        # Predict ADT data
        adt_pred = self.model.predict(X_test)
        # Clip negative values
        adt_pred = np.clip(adt_pred, a_min=0, a_max=None)

        return adt_pred, self.adt_names

    def save(self, path: str, compress: Optional[bool] = True):
        """
        Save the trained pipeline to a file.
        Parameters
        ----------
        path
            Path to the file.
        compress
            Whether to compress the file, default True.
        """
        import joblib
        joblib.dump(self, path, compress=compress)

    def load(self, path: str):
        """
        Load a pretrained pipeline from a file.
        Parameters
        ----------
        path
            Path to the file.
        """
        import joblib
        pretrained_pipe = joblib.load(path)
        self.gex_preprocessor = pretrained_pipe.gex_preprocessor
        self.model = pretrained_pipe.model
        self.gex_names = pretrained_pipe.gex_names
        self.adt_names = pretrained_pipe.adt_names


class ADTPredictorKRREnsemble(ADTPredictor):
    """
    ADT predictor class that uses a kernel ridge regression ensemble model instead of a linear regression model.
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
            Logarithmize data?
            Default 'True' expects non-logarithmized data.
        n_components
            Number of components to use for truncated SVD.
        do_tsvd_before_zscore
            Perform truncated SVD before Z-score normalization?
            Default 'True' works better for downstream training and prediction on data from a single dataset.
            Set to 'False' to extract more robust features that work well across datasets.

        """
        super().__init__(do_log1p, n_components, do_tsvd_before_zscore)
        self.model = KernelRidgeEnsemble()


class KernelRidgeEnsemble:
    """
    Kernel ridge regression ensemble model. Winning model of the NeurIPS 2021 Open Problems in Single-Cell Analysis
    challenge for the modality prediction task from GEX to ADT, proposed by Kaiwen Deng.
    Citation: https://proceedings.mlr.press/v176/lance22a.html
    Code adapted from: https://github.com/openproblems-bio/neurips2021_multimodal_topmethods/tree/main/src/
    predict_modality/methods/Guanlab-dengkw
    """
    def __init__(self):
        self.regressors = []

    def fit(self, X: np.ndarray, y: np.ndarray):
        """
        Fit the kernel ridge regression ensemble.
        Parameters
        ----------
        X
            GEX matrix.
        y
            ADT matrix.
        """
        from tqdm.auto import tqdm
        batch2idxs = {'s1d3': [9185, 14668], 's1d2': [4721, 9184], 's3d7': [55813, 66174], 's3d1': [37254, 45835], 's1d1': [0, 4720],
         's2d1': [14669, 24021], 's2d4': [24022, 29047], 's3d6': [45836, 55812], 's2d5': [29048, 37253]}
        batches_splits = [['s3d7', 's2d4', 's2d1', 's1d2'], ['s1d1', 's2d5', 's1d3', 's3d6', 's3d1'], ['s2d5', 's2d1', 's3d6', 's3d7'],
         ['s1d3', 's1d2', 's1d1', 's2d4', 's3d1'], ['s1d2', 's2d4', 's3d7', 's1d1'],
         ['s2d5', 's2d1', 's3d6', 's3d1', 's1d3'], ['s1d2', 's2d4', 's1d1', 's3d1'],
         ['s3d6', 's1d3', 's3d7', 's2d1', 's2d5'], ['s2d4', 's1d1', 's2d5', 's1d2'],
         ['s3d7', 's2d1', 's3d1', 's3d6', 's1d3']]
        # Fit the model
        for split in tqdm(batches_splits):
            print("initializing")
            # kernel = RBF(length_scale=10)
            length_scale = 10
            gamma = 1 / (2 * length_scale ** 2)
            regressor = KernelRidge(alpha=0.2, kernel='rbf', gamma=gamma)
            mask = np.zeros(X.shape[0], dtype=bool)
            for batch in split:
                mask[batch2idxs[batch][0]:batch2idxs[batch][1]] = True
            print("fitting")
            regressor.fit(X[mask, :], y[mask, :])
            self.regressors.append(regressor)

    def predict(self, X: np.ndarray):
        """
        Predict ADT data.
        Parameters
        ----------
        X
            GEX matrix.
        Returns
        -------
        ADT matrix.
        """
        adt_pred = np.zeros((1000, 134), dtype=np.float32)
        for regressor in self.regressors:
            adt_pred += regressor.predict(X)
        adt_pred /= len(self.regressors)
        np.clip(adt_pred, 0, None, out=adt_pred)
        return adt_pred

    def save(self, path: str, compress: Optional[bool] = True):
        """
        Save the trained model to a file.
        Parameters
        ----------
        path
            Path to the file.
        compress
            Whether to compress the file, default True.
        """
        import joblib
        joblib.dump(self, path, compress=compress)

    def load(self, path: str):
        """
        Load a pretrained model from a file.
        Parameters
        ----------
        path
            Path to the file.
        """
        import joblib
        pretrained_regressors = joblib.load(path)
        self.regressors = pretrained_regressors.regressors

