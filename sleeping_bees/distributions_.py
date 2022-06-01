
import nolds
import pickle
import numpy as np
import pandas as pd
import warnings
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from scipy.stats import skew, kurtosis
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score
warnings.filterwarnings("ignore")

# Get the spatial mean within the square of coordinates x1,x2,y1,y2

# the dimensions are ordered (t, y, x)
# the output are spatially averaged time segments


def squareROI(awake: np.ndarray, sleep: np.ndarray, x1: int, x2: int, y1: int, y2: int, LENGTH: int = 100) -> pd.DataFrame:

    g_awake = []
    for n in range(0, len(awake), LENGTH):
        g_awake.append(awake[n:n+LENGTH, :, :])

    g_sleep = []
    for n in range(0, len(sleep), LENGTH):
        g_sleep.append(sleep[n:n+LENGTH, :, :])

    means = []
    for g_ in g_awake:
        means.append(np.mean(g_[:, y1:y2, x1:x2], axis=(1, 2)))

    for g_ in g_sleep:
        means.append(np.mean(g_[:, y1:y2, x1:x2], axis=(1, 2)))

    y = np.asarray([1] * len(g_awake) + [0] * len(g_sleep), dtype=np.float32)
    # Create features for class 1: awake and class 0: sleep

    X = []
    # sample entropy (sampen):Measures the complexity of a time-series, based on approximate entropy
    # Hurst exponent (hurst_rs)
    # The hurst exponent is a measure of the “long-term memory” of a time series.
    # It can be used to determine whether the time series is more, less, or equally likely
    # to increase if it has increased in previous steps.

    # detrended fluctuation analysis (DFA) (dfa)
    # DFA measures the Hurst parameter H, which is very similar to the Hurst exponent.
    # The main difference is that DFA can be used for non-stationary processes (whose mean and/or variance change over time).
    for m, l in zip(means, y):
        X1 = np.std(m)
        X2 = skew(m)
        X3 = kurtosis(m)
        X4 = nolds.sampen(m)
        X5 = nolds.hurst_rs(m)
        X6 = nolds.dfa(m)
        X.append({'STD': X1, 'Skewness': X2, 'Kurtosis': X3,
                  'Entropy': X4, 'Hurst': X5, 'DFA': X6, 'Class': l})
    dataset = pd.DataFrame(X)

    return dataset


def pair_plot(df: pd.DataFrame, file_name):
    sns.set(style="ticks")
    sns.set_palette("deep")
    sns.pairplot(df.iloc[:, :], hue="Class")
    plt.savefig(file_name)
    plt.show()
    plt.clf()


def build_train_rf_model(X_train: np.ndarray, y_train: np.ndarray, save_to: str):
    rfc = RandomForestClassifier()
    rfc.fit(X_train, y_train)
    with open(save_to, "wb") as f:
        pickle.dump(rfc, f)
    return rfc


def plot_feature_importance(rfc: RandomForestClassifier, cols: list, file_name):
    feat_importances = pd.DataFrame(
        rfc.feature_importances_, index=cols, columns=["Importance"])
    feat_importances.sort_values(
        by='Importance', ascending=False, inplace=True)
    feat_importances.plot(kind='bar', figsize=(8, 6))
    # plt.bar(feat_importances.index.to_numpy(), feat_importances["Importance"].to_numpy())
    plt.savefig(file_name)
    plt.show()
    plt.clf()


def run_distributions(awake: np.ndarray, sleep: np.ndarray, x1=70, x2=80, y1=70, y2=80,
                      LENGTH=100, model_fname="rfc.pkl", pair_plot_fname='pairplot.jpg',
                      feature_importance_fname="rfc_feature_importance.jpg",
                      roc_curve_fname="rfc_roc_curve.jpg"):
    df = squareROI(awake, sleep, x1, x2, y1, y2, LENGTH=LENGTH)
    df.replace(to_replace=np.inf, value=1e6, inplace=True)
    df.dropna(inplace=True)
    pair_plot(df, pair_plot_fname)
    X_train, X_test, y_train, y_test = train_test_split(
        df.drop(columns="Class", inplace=False).to_numpy(), df["Class"].to_numpy(), test_size=0.2, stratify=df["Class"].to_numpy())
    rfc = build_train_rf_model(X_train, y_train, model_fname)
    plot_feature_importance(rfc, df.columns[:-1], feature_importance_fname)
    y_pred = rfc.predict(X_test)
    fpr, tpr, thresholds = roc_curve(y_test, y_pred)
    auc = roc_auc_score(y_test, y_pred)
    
    plt.plot(fpr, tpr, 'o')
    plt.xlabel("false positive rates (Specifity)")
    plt.ylabel("true positive reates (Sensivity)")
    plt.title(f'auc score = {auc}')
    plt.grid(True)
    plt.savefig(roc_curve_fname)
    plt.show()
    plt.clf()
    return df, rfc
