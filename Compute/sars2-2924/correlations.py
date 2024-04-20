import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import rcParams
from seismicrna.core.mu import calc_pearson
from seismicrna.table.load import load_pos_table

# Make editable text in PDF files.
rcParams["pdf.fonttype"] = 42

LABELS = [0, 4, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12]


def calculate(ref_table_file: str, *other_table_file: str):
    correlations = dict()
    ref_table = load_pos_table(Path(ref_table_file))
    for other_table in map(load_pos_table, map(Path, other_table_file)):
        correlation = calc_pearson(
            ref_table.fetch_ratio(rel="Mutated", squeeze=True),
            other_table.fetch_ratio(rel="Mutated", squeeze=True)
        )
        sample = other_table.sample
        name, _ = sample.split("-")
        if name == "ctrl2":
            label = 0
        elif name.startswith("aso"):
            label = int(name[len("aso"):])
        else:
            raise ValueError(f"Invalid name: {repr(name)}")
        correlations[label] = correlation
    return pd.Series(correlations).sort_index()


def graph(graph_file: str, correlations: pd.Series):
    fig, (axt, axb) = plt.subplots(nrows=2, sharex="all")
    x = np.arange(len(LABELS))
    for ax in [axt, axb]:
        ax.bar(x, correlations.loc[LABELS])
        ax.spines[["top", "right", "bottom"]].set_visible(False)
    axb.set_xticks(x)
    axb.set_xticklabels(LABELS)
    # Break the y-axis and remove the area between 0.60 and 0.78
    axt.set_ylim(0.78, 1.00)
    axb.set_ylim(0.00, 0.60)
    axt.set_aspect(20.)
    axb.set_aspect(3.)
    fig.subplots_adjust(hspace=0.0)
    plt.savefig(graph_file)


def main(out_file: Path, ref_table_file: Path, *other_table_file: Path):
    graph_file = out_file.with_suffix(".pdf")
    data_file = out_file.with_suffix(".csv")
    data = calculate(ref_table_file, *other_table_file)
    data.to_csv(data_file, header=False)
    graph(graph_file, data)


if __name__ == "__main__":
    main(*map(Path, sys.argv[1:]))
