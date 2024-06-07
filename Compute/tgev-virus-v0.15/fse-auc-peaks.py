import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


TOPLINE = 0.95
SLIPPERY = 12338  # TGEV
#SLIPPERY = 13483  # SARS2
LDI = 13200  # TGEV
#LDI = 14665  # SARS-2
UPSTREAM = 0
BINS = 101
COLOR_SLIPPERY = "#ff0000"
COLOR_LDI = "#0000ff"
COLOR_OTHER = "#e5e5e5"


def read_auc(auc_file: str):
    auc = pd.read_csv(auc_file, header=[0, 1], index_col=[0, 1]).squeeze()
    assert isinstance(auc, pd.Series)
    return auc


def call_tops_dips(auc: pd.Series):
    auc = auc.loc[~np.isnan(auc)]
    positions = auc.index.get_level_values("Position")
    is_top = auc.values >= TOPLINE
    is_dip = auc.values < TOPLINE
    top5 = is_top & np.hstack([[True], is_dip[:-1]])
    top3 = is_top & np.hstack([is_dip[1:], [True]])
    dip5 = is_dip & np.hstack([[True], is_top[:-1]])
    dip3 = is_dip & np.hstack([is_top[1:], [True]])
    tops = list(zip(positions[top5], positions[top3], strict=True))
    dips = list(zip(positions[dip5], positions[dip3], strict=True))
    return tops, dips


def get_widths(index: pd.Index):
    end5s = index.get_level_values(0)
    end3s = index.get_level_values(1)
    return pd.Series(end3s - end5s + 1, index)


def get_color(ends: tuple[int, int]):
    end5, end3 = ends
    if end5 <= SLIPPERY <= end3:
        return COLOR_SLIPPERY
    elif end5 <= LDI <= end3:
        return COLOR_LDI
    return COLOR_OTHER


def get_colors(index: pd.Index):
    return pd.Series(list(map(get_color, index)), index)


def compute_dip_areas(auc: pd.Series):
    tops, dips = call_tops_dips(auc)
    areas = pd.Series({
        (end5, end3): np.nansum(np.minimum(auc.loc[end5: end3] - TOPLINE, 0.))
        for end5, end3 in tops + dips
    })
    return areas.loc[areas.index[np.argsort(areas)]]


def graph_dip_hist(areas: pd.Series, graph_file):
    fig, ax = plt.subplots()
    ax.hist(-areas, bins=BINS)
    ax.set_xlabel("Area of Dip")
    ax.set_ylabel("Number of Dips")
    ax.set_yscale("log")
    plt.savefig(graph_file)
    plt.close()


def graph_dip_ranks(areas: pd.Series, graph_file):
    fig, ax = plt.subplots()
    ax.scatter(np.arange(areas.size), areas, c=get_colors(areas.index))
    ax.set_xlabel("Rank by Area")
    ax.set_ylabel("Area of Dip")
    plt.savefig(graph_file)
    plt.close()


def graph_dip_widths(index: pd.Index, graph_file):
    widths = get_widths(index)
    fig, ax = plt.subplots()
    ax.hist(widths, bins=BINS)
    ax.set_xlabel("Width of Dip")
    ax.set_ylabel("Number of Dips")
    plt.savefig(graph_file)
    plt.close()


def main(auc_file: str, output_prefix: str):
    output_prefix = Path(output_prefix)
    areas = compute_dip_areas(read_auc(auc_file))
    areas.to_csv(output_prefix.with_suffix(".csv"), header=False)
    graph_dip_ranks(areas, output_prefix.with_suffix(".pdf"))
    graph_dip_widths(areas.index, output_prefix.with_suffix(".pdf"))


if __name__ == "__main__":
    main(*sys.argv[1:])
