import gc

import anndata as ad
from evaluate import evaluate
from prediction import ADTPredictorKRREnsemble, ADTPredictor
import pandas as pd
import os
import random
import numpy as np
from tqdm.auto import tqdm
import argparse

# parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("-d", "--datasets_path", help="datasets path containing the directory "
                                                  "openproblems_bmmc_cite_phase2_rna",
                    default="../../../../../PycharmProjects/ModalityPrediction/datasets")
parser.add_argument("-cv", "--cross_validation", help="run cross validation instead of training on all data",
                    action="store_true")
parser.add_argument("-lr", "--linear_regression", help="use linear regression instead of "
                                                       "kernel ridge regression ensemble", action="store_true")
args = parser.parse_args()

# seed everything
seed = 0
random.seed(seed)
np.random.seed(seed)
os.environ['PYTHONHASHSEED'] = str(seed)

# load data
dataset_path = args.datasets_path + "/openproblems_bmmc_cite_phase2_rna/openproblems_bmmc_cite_phase2_rna" \
                                    ".censor_dataset.output_"
gex_train = ad.read_h5ad(dataset_path + "train_mod1.h5ad")
gex_test = ad.read_h5ad(dataset_path + "test_mod1.h5ad")
adt_train = ad.read_h5ad(dataset_path + "train_mod2.h5ad")
adt_test = ad.read_h5ad(dataset_path + "test_mod2.h5ad")

def batch2idxs(gex):
    batch2idxs = dict()
    for batch in gex.obs["batch"].cat.categories.tolist():
        all_idxs = gex.obs.batch.isin([batch]).to_numpy().nonzero()[0]
        batch2idxs[batch] = (all_idxs.min(), all_idxs.max() + 1)
    return batch2idxs

if args.cross_validation:
    # run out-of-site cross validation
    gex_all = ad.concat((gex_train, gex_test), join='outer')
    adt_all = ad.concat((adt_train, adt_test), join='outer')
    batch_labels = gex_all.obs["batch"].cat.categories.tolist()     # s1d1, s1d2, s2d1, ..
    site_labels = list(set([label[1] for label in batch_labels]))   # 1, 2, ..
    cv_scores = []
    for i, site in tqdm(enumerate(site_labels)):
        # train on all data except from site
        train_batch_labels = [label for label in batch_labels if label[1] != site]
        valid_batch_labels = [label for label in batch_labels if label[1] == site]
        print(f"Training on {train_batch_labels}, validating on {valid_batch_labels}")
        train_indices = gex_all.obs["batch"].isin(train_batch_labels)
        valid_indices = gex_all.obs["batch"].isin(valid_batch_labels)
        X_train = gex_all[train_indices].X.toarray()
        X_valid = gex_all[valid_indices].X.toarray()
        Y_train = adt_all[train_indices].X.toarray()
        Y_valid = adt_all[valid_indices].X.toarray()
        if args.linear_regression:
            pipe = ADTPredictor(do_log1p=False)
        else:
            pipe = ADTPredictorKRREnsemble(do_log1p=False, batch2idxs=batch2idxs(gex_all[train_indices]))
        pipe.fit(gex_train=X_train,
                 adt_train=Y_train,
                 gex_test=X_valid)
        Y_pred, adt_names = pipe.predict(X_valid)
        rmse, pearson, spearman = evaluate(Y_pred, Y_valid)
        cv_scores.append((site, rmse, pearson, spearman))
        del X_train, X_valid, Y_train, Y_valid, Y_pred, pipe
        gc.collect()

    cv_scores_df = pd.DataFrame(data=cv_scores, columns=["site", "rmse", "pearson", "spearman"])
    print(cv_scores_df)
    for metric in ["rmse", "pearson", "spearman"]:
        print(cv_scores_df.agg({metric: ["mean", "std"]}))
else:
    # train on all data and save the pipeline
    if args.linear_regression:
        pipe = ADTPredictor(do_log1p=False)
        pipeline_name = "ADTPredictor"
    else:
        pipe = ADTPredictorKRREnsemble(do_log1p=False, batch2idxs=batch2idxs(gex_train))
        pipeline_name = "ADTPredictorKRREnsemble"
    pipe.fit(gex_train=gex_train.X.toarray(), adt_train=adt_train.X.toarray(), gex_test=gex_test.X.toarray(),
             gex_names=gex_train.var_names.to_numpy(), adt_names=adt_train.var_names.to_numpy())
    pipe.save(f"{pipeline_name}_neuripstrain_alltypes.joblib")
    adt_pred, adt_names = pipe.predict(gex_test.X.toarray())
    evaluate(adt_pred, adt_test.X.toarray())


