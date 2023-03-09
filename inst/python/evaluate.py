from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import mean_squared_error


def evaluate(y_pred, y_test, verbose=True):
    rmse = mean_squared_error(y_pred, y_test, squared=False)

    pearson_sum = 0
    spearman_sum = 0
    for i in range(len(y_test)):
        pearson_sum += pearsonr(y_test[i], y_pred[i])[0]
        spearman_sum += spearmanr(y_test[i], y_pred[i])[0]

    pearson_corr = pearson_sum / len(y_test)
    spearman_corr = spearman_sum / len(y_test)
    if verbose:
        print("RMSE:", rmse)
        print("Pearson correlation:", pearson_corr)
        print("Spearman correlation:", spearman_corr)
    return rmse, pearson_corr, spearman_corr
