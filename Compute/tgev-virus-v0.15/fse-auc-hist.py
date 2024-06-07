import sys

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


SLIPPERY = 12338
NUP = 100
NDN = 0
BINS = np.linspace(0., 1., 101)
WIDTH = 0.004


def read_auc(auc_file: str):
    auc = pd.read_csv(auc_file, header=[0, 1], index_col=[0, 1]).squeeze()
    assert isinstance(auc, pd.Series)
    return auc


def compute_histograms(bins: np.ndarray, auc: pd.Series, nup: int, ndn: int):
    end5 = SLIPPERY - nup
    end3 = SLIPPERY + ndn
    auc_in = auc.loc[end5: end3]
    auc_out = pd.concat([auc.loc[: end5 - 1], auc.loc[end3 + 1:]])
    hist_in, _ = np.histogram(auc_in, bins=bins)
    hist_out, _ = np.histogram(auc_out, bins=bins)
    return hist_in, hist_out


def graph_histograms(bins: np.ndarray,
                     hist_in: np.ndarray,
                     hist_out: np.ndarray,
                     graph_file: str):
    x = (bins[:-1] + bins[1:]) / 2.
    fig, ax_out = plt.subplots()
    color_out = "#0072b2"
    ax_out.set_xlabel("AUC-ROC")
    ax_out.set_ylabel("Number of Positions", color=color_out)
    ax_out.bar(x, hist_out, width=-WIDTH, align="edge", label="Away from FSE", color=color_out)
    ax_out.tick_params(axis="y", labelcolor=color_out)
    ax_in = ax_out.twinx()
    color_in = "#e69f00"
    ax_in.set_ylabel("Number of Positions", color=color_in)
    ax_in.bar(x, hist_in, width=WIDTH, align="edge", label="Near FSE", color=color_in)
    ax_in.tick_params(axis="y", labelcolor=color_in)
    fig.tight_layout()
    fig.legend(loc=2)
    plt.savefig(graph_file)
    plt.close()


def main(auc_file: str, graph_file: str):
    auc = read_auc(auc_file=auc_file)
    hist_in, hist_out = compute_histograms(BINS, auc, NUP, NDN)
    graph_histograms(BINS, hist_in, hist_out, graph_file)


if __name__ == "__main__":
    main(*sys.argv[1:])
