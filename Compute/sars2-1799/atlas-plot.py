import sys
from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import colormaps
from matplotlib.cm import ScalarMappable
from matplotlib.patches import Rectangle

from seismicrna.core.header import format_clust_name, list_clusts
from seismicrna.core.path import fill_whitespace
from seismicrna.core.rna import RNAProfile, RNAState, RNAStructure, from_ct
from seismicrna.table.base import R_ADJ_TITLE
from seismicrna.table.load import ClustPosTableLoader, ClustFreqTableLoader

CMAP_NAME = "coolwarm"
COLOR_MAP = colormaps.get_cmap(CMAP_NAME)


def compute_auc(structs: Iterable[RNAStructure],
                profiles: Iterable[RNAProfile]):
    """ Compute the AUC-ROC for each profile vs. each structure. """
    auc = defaultdict(dict)
    for struct, profile in product(structs, profiles):
        state = RNAState.from_struct_profile(struct, profile)
        auc[state.data_name][state.title] = state.auc
    return pd.DataFrame.from_dict(auc)


def compute_proportion(freq_table: ClustFreqTableLoader):
    """ Compute the cluster proportions. """
    proportions = dict()
    for order in freq_table.header.orders:
        freqs = freq_table.data.loc[R_ADJ_TITLE, order]
        props = freqs / freqs.sum()
        for clust in list_clusts(order):
            name = fill_whitespace(format_clust_name(order, clust), fill="-")
            proportions[name] = props.loc[clust]
    return pd.Series(proportions)


def sort_clusters(auc: pd.DataFrame, prop: pd.Series):
    """ Sort clusters in descending order by their proportions. """
    prop_use = prop.loc[auc.columns]
    prop_ordered = prop_use.sort_values(ascending=False)
    auc_ordered = auc.loc[:, prop_ordered.index]
    return auc_ordered, prop_ordered


def get_color(auc: float):
    """ Get the color for an AUC value. """
    return COLOR_MAP(1. - auc)


def graph_auc(graph_file: Path, auc: pd.DataFrame, prop: pd.Series):
    num_structs, num_clusts = auc.shape
    # Compute the coordinates.
    values, widths = sort_clusters(auc, prop)
    if widths.min() < 0.:
        raise ValueError(f"Expected proportions to be ≥ 0, but got {widths}")
    rights = np.cumsum(widths)
    if not np.isclose(total := rights.iloc[-1], 1.):
        raise ValueError(f"Expected proportions to sum to 1, but got {total}")
    lefts = rights - widths
    heights = pd.Series(1., index=auc.index[::-1])
    tops = np.cumsum(heights)
    bottoms = tops - heights
    # Set up the axis.
    fig, ax = plt.subplots()
    fig.colorbar(ScalarMappable(cmap=COLOR_MAP), ax=ax)
    ax.set_xlim(lefts.iloc[0], rights.iloc[-1])
    ax.set_ylim(bottoms.iloc[0], tops.iloc[-1])
    xticks = (rights + lefts) / 2.
    yticks = (bottoms + tops) / 2.
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks.index)
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks.index)
    # Plot the values.
    for clust, clust_values in values.items():
        for struct, value in clust_values.items():
            ax.add_patch(Rectangle(
                (lefts.loc[clust], bottoms.loc[struct]),
                widths.loc[clust],
                heights.loc[struct],
                facecolor=get_color(value)
            ))
    plt.savefig(graph_file)
    return values


def main(out_file: Path,
         ct_file: Path,
         pos_table_file: Path,
         freq_table_file: Path,
         order: int):
    graph_file = out_file.with_suffix(".pdf")
    csv_file = out_file.with_suffix(".csv")
    values = graph_auc(
        graph_file,
        compute_auc(
            from_ct(ct_file),
            ClustPosTableLoader(pos_table_file).iter_profiles(order=order)
        ),
        compute_proportion(ClustFreqTableLoader(freq_table_file))
    )
    values.to_csv(csv_file)


if __name__ == "__main__":
    out_file, ct_file, pos_table_file, freq_table_file, order = sys.argv[1:]
    main(Path(out_file),
         Path(ct_file),
         Path(pos_table_file),
         Path(freq_table_file),
         int(order))
