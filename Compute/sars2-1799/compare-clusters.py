from collections import defaultdict
from itertools import product

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import pearsonr


def make_heatmap(aso: int, order: int):
    file_in = f"out/LNA{aso}-rep1__and__LNA{aso}-rep2/graph/sars2/fse/scatter_clustered-{order}-x_m-ratio-q0.csv"
    data = pd.read_csv(file_in, header=list(range(4)), index_col=list(range(2)))
    rep1 = data[f"LNA{aso}-rep1"]["Mutated"][str(order)]
    rep2 = data[f"LNA{aso}-rep2"]["Mutated"][str(order)]
    corr = defaultdict(dict)
    for cluster1, cluster2 in product(range(1, order + 1), repeat=2):
        series1 = rep1[str(cluster1)]
        series2 = rep2[str(cluster2)]
        mask = np.isnan(series1) | np.isnan(series2)
        corr[cluster1][cluster2] = pearsonr(series1[~mask], series2[~mask]).statistic
    corr = pd.DataFrame(corr)
    csv_file_out = f"compare-clusters/LNA{aso}_{order}.csv"
    corr.to_csv(csv_file_out)
    pdf_file_out = f"compare-clusters/LNA{aso}_{order}.pdf"
    plt.imshow(corr, cmap="viridis", vmin=0., vmax=1., aspect=1.)
    if order == 1:
        plt.colorbar()
    plt.savefig(pdf_file_out)
    plt.close()


def make_heatmaps():
    for aso in range(10):
        for order in range(10):
            try:
                make_heatmap(aso, order)
                print(aso, order)
            except FileNotFoundError:
                pass


if __name__ == "__main__":
    make_heatmaps()

