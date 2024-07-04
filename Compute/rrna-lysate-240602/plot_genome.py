import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Arc

from seismicrna.core.rna import RNAStructure, from_ct
from seismicrna.table.load import load_pos_table


def mm_to_inch(mm: float):
    """ Convert millimeters to inches. """
    return mm / 25.4


def inch_to_point(inch: float):
    """ Convert inches to points. """
    return inch * 72.


def mm_to_point(mm: float):
    """ Convert millimeters to points. """
    return inch_to_point(mm_to_inch(mm))
    

plt.rcParams["font.family"] = "Helvetica Neue"
plt.rcParams["font.size"] = 6.
plt.rcParams["svg.fonttype"] = "none"
AXIS_ASPECT = 500.
ARC_ASPECT = 1.


def round_inc(x: float | int, inc: float | int, method: str = "nearest"):
    if inc == 0:
        return x
    if method == "nearest":
        return round(x / inc) * inc
    if method == "floor":
        return math.floor(x / inc) * inc
    if method == "ceil":
        return math.ceil(x / inc) * inc
    raise ValueError(f"Invalid method: {repr(method)}")


def load_structure(ct_file: Path):
    """ Load the structure from a CT file. """
    structures = list(from_ct(ct_file))
    if len(structures) != 1:
        raise ValueError(
                f"Expected exactly 1 structure, but got {len(structures)}"
        )
    return structures[0]


def compute_depths(structure: RNAStructure):
    """ Compute the depth of each base in the structure. """
    positions = structure.section.range_int
    opens = np.logical_and(structure.is_paired, structure.table > positions)
    closes = np.logical_and(structure.is_paired, structure.table < positions)
    depths = np.cumsum(np.asarray(opens, dtype=int)
                       - np.hstack([np.array([0]),
                                    np.asarray(closes[:-1], dtype=int)]))
    assert np.all(depths >= 0)
    return pd.Series(depths, positions)


def partition_structure(structure: RNAStructure, nrows: int):
    """ Partition an RNA structure into rows. """
    end5 = structure.section.end5
    end3 = structure.section.end3
    # Target locations of row breaks.
    targets = np.asarray(np.round(np.linspace(end5, end3 + 1, nrows + 1)[1: -1]),
                         dtype=int)
    assert len(targets) == nrows - 1
    # Locations of gap between stems nearest each row break.
    depths = compute_depths(structure)
    gaps = depths.index.values[depths == 0]
    nearest_gaps = list()
    for target in targets:
        distances = np.abs(gaps - target)
        nearest_gaps.append(int(gaps[np.min(distances) == distances][0]))
    # Range for each row.
    if nrows > 1:
        partitions = [(end5, nearest_gaps[0])]
        for gap in range(nrows - 2):
            partitions.append((nearest_gaps[gap] + 1, nearest_gaps[gap + 1]))
        partitions.append((nearest_gaps[nrows - 2] + 1, end3))
    else:
        partitions = [(end5, end3)]
    return partitions


def plot_rows(structure: RNAStructure,
              trend: pd.Series,
              mutated: pd.Series | None,
              output_file: Path,
              nrows: int = 1):
    """ Plot the structure, mutation rate, and AUC-ROC for each row. """
    if mutated is not None:
        mutated = mutated.reindex(trend.index)
    partitions = partition_structure(structure, nrows)
    fig, axes = plt.subplots(nrows=len(partitions),
                                       sharex=False,
                                       sharey=True)
    fig.set_size_inches(mm_to_inch(180.),
                        mm_to_inch(240.))
    fig.set_layout_engine("compressed")
    axes = np.atleast_1d(axes)
    axes[-1].set_xlabel("Position")
    axes[-1].set_ylabel("Attribute")
    for ax, (row5, row3) in zip(axes, partitions, strict=True):
        print(f"Row: {row5} - {row3}")
        positions = structure.section.range[row5: row3].get_level_values("Position")
        # Axis setup.
        ax.set_xlim(row5 - 0.5, row3 + 0.5)
        ax.set_xticks(np.arange(round_inc(row5, 500, method="ceil"),
                                round_inc(row3, 500, method="floor") + 1,
                                500),
                      minor=False)
        ax.set_xticks(np.arange(round_inc(row5, 100, method="ceil"),
                                round_inc(row3, 100, method="floor") + 1,
                                100),
                      minor=True)
        ax.set_ylim(0., 1.)
        ax.set_yticks(np.linspace(0., 1., 3), minor=False)
        ax.set_yticks(np.linspace(0., 1., 11), minor=True)
        ax.set_yticklabels([], minor=True)
        ax.tick_params(axis="both", which="both", length=0.)
        ax.set_aspect(AXIS_ASPECT)
        ax.spines.top.set_visible(False)
        ax.spines.right.set_visible(False)
        ax.spines.left.set_linewidth(mm_to_point(0.2))
        ax.spines.bottom.set_linewidth(mm_to_point(0.2))
        ax.grid(axis="x",
                which="both",
                visible=True,
                linewidth=mm_to_point(0.2),
                color="#f2f2f2")
        ax.grid(axis="y",
                which="both",
                visible=True,
                linewidth=mm_to_point(0.2),
                color="#e5e5e5")
        # Mutation rate.
        if mutated is not None:
            ax.bar(positions,
                   mutated.loc[positions],
                   width=1.,
                   color="#d55e00")
        # Trend.
        ax.plot(positions,
                trend.loc[positions],
                linewidth=mm_to_point(0.2),
                color="#009e73")
        ax.fill_between(positions,
                        trend.loc[positions],
                        y2=1.,
                        color="#83e8cc")
        # Secondary structure.
        for position in positions:
            partner = int(structure.table.at[position].iloc[0])
            if partner > position:
                if not row5 <= partner <= row3:
                    raise ValueError(f"Partner {partner} of base {position} "
                                     f"is out of bounds for row {row5} - {row3}")
                x = (position + partner) / 2.
                width = partner - position
                height = width * ARC_ASPECT / AXIS_ASPECT
                ax.add_patch(Arc((x, 0.),
                                 width,
                                 height,
                                 theta1=0.,
                                 theta2=180.,
                                 linewidth=mm_to_point(0.03),
                                 edgecolor="#56b4e9"))
    plt.savefig(output_file)
    plt.close()


def run(ct_file: str,
        trend_file: str,
        table_file: str | None,
        output_file: str,
        nrows: int = 1):
    structure = load_structure(Path(ct_file))
    trend = pd.read_csv(Path(trend_file),
                        index_col=[0, 1],
                        header=[0, 1])
    if trend.columns.size > 1:
        raise ValueError(f"Expected 1 column, but got {trend.columns.size}")
    trend = trend.iloc[:, 0]
    assert isinstance(trend, pd.Series)
    if table_file:
        table = load_pos_table(Path(table_file))
        mutated = table.fetch_ratio(rel="Mutated", squeeze=True)
    else:
        mutated = None
    plot_rows(structure,
              trend,
              mutated,
              Path(output_file),
              nrows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("ct", type=str)
    parser.add_argument("trend", type=str)
    parser.add_argument("output", type=str)
    parser.add_argument("--table", "-t", type=str)
    parser.add_argument("--nrows", "-n", type=int, default=1)
    args = parser.parse_args()
    run(args.ct,
        args.trend,
        args.table,
        args.output,
        args.nrows)
