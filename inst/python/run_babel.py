import anndata as ad

from evaluate import evaluate
from prediction import ADTPredictorBabel

dataset_path = "../../../../../PycharmProjects/ModalityPrediction/datasets/openproblems_bmmc_cite_phase2_rna/openproblems_bmmc_cite_phase2_rna.censor_dataset.output_"
gex_train = ad.read_h5ad(dataset_path + "train_mod1.h5ad")
gex_test = ad.read_h5ad(dataset_path + "test_mod1.h5ad")
adt_train = ad.read_h5ad(dataset_path + "train_mod2.h5ad")
adt_test = ad.read_h5ad(dataset_path + "test_mod2.h5ad")

# using high-level interface
# set use_vanilla_nn to True to use a Vanilla NN instead of BabelDance
pipe = ADTPredictorBabel(do_log1p=False, use_vanilla_nn=False)
# fit on training data
# gex_test is optional and is used for transductive preprocessing if provided
# gex_names and adt_names are optional and should refer to the variable names of gex_train and adt_train
# if not provided, the predict() method will assume that all the columns of the test GEX matrix are in the same order as the training GEX matrix

pipe.fit(gex_train=gex_train.X.toarray(), adt_train=adt_train.X.toarray(), gex_test=gex_test.X.toarray(),
         gex_names=gex_train.var_names.to_numpy(), adt_names=adt_train.var_names.to_numpy())

adt_pred, adt_names = pipe.predict(gex_test.X.toarray())
evaluate(adt_pred, adt_test.X.toarray())

# adt names are also stored as a property
print(pipe.adt_names)

# save the trained pipeline to a file
pipe.save("ADTPredictorBabel_neuripstrain_alltypes.joblib")