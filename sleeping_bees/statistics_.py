import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import warnings
# from statannot import add_stat_annotation

sns.set(color_codes=True)
warnings.filterwarnings("ignore")


def plot_feature_distribution(data: pd.DataFrame, save_file_fmt="feature_distributions_{}.jpg"):
    a = data.Class.unique()
    my_pal = {Class: "yellow" if Class ==
              a[0] else "lime" for Class in a}
    for c in data.columns[:-1]:
        ax = sns.boxplot(
            x="Class", y=c, data=data, palette=my_pal)
        # add_stat_annotation(ax=ax, data=data, x="Class", y=c, box_pairs=[("Awake", "Sleep")],
        #                     test='t-test_paired', comparisons_correction=None, text_format='star',
        #                     loc='inside', verbose=1)
        plt.savefig(save_file_fmt.format(c))
        plt.show()
        plt.clf()

    plt.rcParams["figure.figsize"] = (20, 20)
    sns.set(style="ticks")
    sns.pairplot(data.iloc[:, :], hue="Class", palette=my_pal)
    plt.savefig(save_file_fmt.format('all'))

    plt.show()
    