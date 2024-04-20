import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Arc
from seismicrna.core.rna import RNAStructure, from_ct
from seismicrna.core.seq import POS_NAME
from seismicrna.mask.report import MaskReport

# Report file containing the section information.
FSE_REPORT_FILE = Path("out/LNA0-pool/mask/sars2/fse/mask-report.json")


# Coordinates of the main row.
MAIN5 = 252
MAIN3 = 423


# Axis setup
YMIN = 0.5
YMAX = 1.25


def get_fse_coords(report_file: Path):
    """ Get the 5' and 3' coordinates of the FSE section. """
    report = MaskReport.load(report_file)
    assert report.sect == "fse"
    return report.end5, report.end3


def get_struct(struct_file: Path):
    """ Get the first structure in a CT file. """
    for struct in from_ct(struct_file):
        return struct
    raise ValueError(f"{struct_file} contained no structures")


def draw_base_pairs(ax,
                    struct: RNAStructure,
                    end5: int,
                    end3: int,
                    lna5: int,
                    lna3: int):
    """ Draw the base pairs in the structure between end5 and end3. """
    assert end5 <= MAIN5 <= MAIN3 <= end3
    # Select the base pairs involving the main row.
    pairs = struct.table.loc[MAIN5: MAIN3]
    pairs.index = pairs.index.get_level_values(POS_NAME)
    # Remove unpaired bases.
    pairs = pairs[pairs != 0]
    # Determine the ranges of the pairs before and after the main row.
    before = pairs[pairs < MAIN5]
    if before.size > 0:
        raise ValueError("Got base pairs extending backwards before main row")
    after = pairs[pairs > MAIN3]
    if after.size == 0:
        raise ValueError("Got no base pairs extending onwards after main row")
    # Map positions in the RNA to coordinates in the graph.
    after_min = after.min()
    after_max = after.max()

    def map_position(position: int | float):
        if end5 <= position <= end3:
            return float(position), 1.
        return float(end5 + (after_max - position)), YMAX

    def plot_pair(p1: int, p2: int):
        if lna5 <= p1 <= lna3 or lna5 <= p2 <= lna3:
            # The LNA overlaps the base pair.
            color = "#ff0000"
        else:
            # The LNA does not overlap the base pair.
            color = "#e5e5e5"
        if end5 <= p1 <= end3 and end5 <= p2 <= end3:
            # The pairs are both on the main row: connect with an arc.
            width = float(abs(p2 - p1))
            height = width * 0.001
            center = map_position((p1 + p2) / 2.)
            ax.add_patch(Arc(center,
                             width=width,
                             height=height,
                             theta1=0.,
                             theta2=180.,
                             color=color))
        else:
            # The pairs are on different rows: connect with a line.
            x1, y1 = map_position(p1)
            x2, y2 = map_position(p2)
            ax.plot([x1, x2], [y1, y2], color=color)

    for p1, p2 in set(pairs.items()):
        plot_pair(p1, p2)

    # Draw the LNA above the structure.
    ax.plot([map_position(lna5), map_position(lna3)], [YMAX, YMAX], color="#00ff00")
    # Set the limits of the axis.
    print("Setting x limits",
          (after_max, map_position(after_max)[0], end5),
          (after_min, map_position(after_min)[0], end3))
    ax.set_xlim(min(map_position(after_max)[0], end5),
                max(map_position(after_min)[0], end3))
    ax.set_ylim(YMIN, YMAX)
    ax.set_yticks(np.linspace(YMIN, 1.0, 6))


def draw_correlation(ax, correl: pd.DataFrame):
    positions = correl.index.get_level_values(POS_NAME)
    values = np.maximum(correl.values, YMIN)
    ax.fill_between(positions, values, 1., color="#aaaaff")
    ax.plot(positions, values, color="#0000ff")


def get_lna_coords(lna_file: str, lna: str):
    lnas_coords = pd.read_csv(lna_file, index_col="LNA")
    lna5 = int(lnas_coords.loc[lna, "Target 5' End (in 1799 nt)"])
    lna3 = int(lnas_coords.loc[lna, "Target 3' End (in 1799 nt)"])
    return lna5, lna3


def main(struct_file: str,
         lna_file: str,
         lna: str,
         correl_file: str,
         graph_file: str):
    fig, ax = plt.subplots()
    correl = pd.read_csv(correl_file, index_col=[0, 1], header=[0, 1]).squeeze()
    assert isinstance(correl, pd.Series)
    struct = get_struct(struct_file)
    end5, end3 = get_fse_coords(FSE_REPORT_FILE)
    lna5, lna3 = get_lna_coords(lna_file, lna)
    draw_correlation(ax, correl)
    draw_base_pairs(ax, struct, end5, end3, lna5, lna3)
    ax.grid(axis="x", visible=True, color="#f2f2f2", linewidth=0.25)
    ax.grid(axis="y", visible=True, color="#e5e5e5", linewidth=0.5)
    ax.spines[["left", "top", "right", "bottom"]].set_visible(False)
    plt.setp(ax.get_xticklines(), visible=False)
    plt.setp(ax.get_yticklines(), visible=False)
    plt.tight_layout()
    plt.savefig(graph_file)


if __name__ == "__main__":
    main(*sys.argv[1:])
