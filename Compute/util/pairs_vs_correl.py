#!/Users/mfa/conda/envs/seismic/bin/python

""" Graph the correlation below the predicted RNA structure. """

from argparse import ArgumentParser
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd
from matplotlib import patches
from matplotlib import pyplot as plt

from seismicrna.core import path, rna, seq
from seismicrna.table.load import load_pos_table

# Make editable text in PDF files.
mpl.rcParams["pdf.fonttype"] = 42

# Set the height-to-width ratio of each arc.
GOLDEN_RATIO = (1. + 5. ** 0.5) / 2.
STRUCT_CORR_RATIO = 0.5

CORR_GRID = 6
CORR_THRESH = 0.9
CORR_MAX = 1.
CORR_MIN = 0.

POSITION_STEP = 300

PAIR_NONE = "#e5e5e5"
PAIR_TRUE = "#0072b2"
PAIR_FALSE = "#e5e5e5"
ASO_COLOR = "#f0e442"
NO_ASO_COLOR = "#009e73"


def calc_arc_aspect_ratio(ax: plt.Axes):
    """ Calculate the width-to-height ratio that will produce the proper
    ratio of height of the correlation to height of the structure. """
    xmin, xmax = ax.get_xlim()
    xsize = xmax - xmin
    ysize = xsize / 2.
    return ysize / STRUCT_CORR_RATIO


def parse_structs_path(structs_path: Path):
    return path.parse(structs_path, *path.CT_FILE_SEGS)


def parse_correl_path(correl_path: Path):
    return path.parse(correl_path, *path.SECT_DIR_SEGS, path.GraphSeg)


def get_ref_name(structs_path: Path, correl_path: Path):
    structs_fields = parse_structs_path(structs_path)
    structs_ref = str(structs_fields[path.REF])
    correl_fields = parse_correl_path(correl_path)
    correl_ref = str(correl_fields[path.REF])
    if structs_ref != correl_ref:
        raise ValueError("References from structure and CORR differ: "
                         f"{structs_ref} ≠ {correl_ref}")
    return structs_ref


def get_table_path_from_structs_path(structs_path: Path):
    structs_fields = parse_structs_path(structs_path)
    table_fields = dict()
    table_fields[path.TOP] = structs_fields[path.TOP]
    table_fields[path.SAMP] = structs_fields[path.SAMP]
    table_fields[path.CMD] = path.TABLE
    table_fields[path.REF] = structs_fields[path.REF]
    table_fields[path.SECT] = structs_fields[path.PROFILE].split("__")[0]
    table_fields[path.TABLE] = path.MASK_TABLE
    table_fields[path.EXT] = path.CSV_EXT
    return path.build(*path.POS_TABLE_SEGS, **table_fields)


def iter_rna_structs(structs_path: Path):
    """ Yield the RNA structures. """
    # Make an RNA state from each structure.
    for struct in rna.from_ct(structs_path):
        yield struct


def iter_rna_states(structs_path: Path):
    """ Yield the RNA states. """
    # Load the RNA profile (mutation rates).
    table = load_pos_table(get_table_path_from_structs_path(structs_path))
    profiles = list(table.iter_profiles())
    if len(profiles) != 1:
        raise ValueError(f"Expected 1 profile, but got {profiles}")
    profile = profiles[0]
    # Make an RNA state from each structure.
    for struct in rna.from_ct(structs_path):
        yield rna.RNAState.from_struct_profile(struct, profile)


def load_correl_data(correl_path: Path):
    """ Load the normalized RMSD data from a file. """
    return pd.read_csv(correl_path, header=[0, 1], index_col=[0, 1])


def get_aso_coords(asos_path: Path, ref: str):
    """ Get the coordinates of an ASO. """
    coords = pd.read_csv(asos_path, index_col="Reference")
    end5 = int(coords.loc[ref, "5' End"])
    end3 = int(coords.loc[ref, "3' End"])
    return end5, end3


def color_base_pair(pair: tuple[int, int],
                    aso_coords: tuple[int, int],
                    correl: pd.Series):
    """ Determine the color of a base pair. """
    # Determine whether the 5' and 3' bases in the pair are targeted by
    # the ASO.
    p5, p3 = pair
    a5, a3 = aso_coords
    target5 = a5 <= p5 <= a3
    target3 = a5 <= p3 <= a3
    # Choose a color based on which bases are targeted and the CORR.
    if target5:
        if target3:
            # The pair lies entirely within the ASO-targeted region.
            return PAIR_NONE
        # The ASO targets only the 5' base.
        d = float(correl.at[p3])
        if d <= CORR_THRESH:
            return PAIR_TRUE
        if d > CORR_THRESH:
            return PAIR_FALSE
        return PAIR_NONE
    if target3:
        # The ASO targets only the 3' base.
        d = float(correl.at[p5])
        if d <= CORR_THRESH:
            return PAIR_TRUE
        if d > CORR_THRESH:
            return PAIR_FALSE
        return PAIR_NONE
    # The pair does not involve a base targeted by the ASO.
    return PAIR_NONE


def graph_structure(ax: plt.Axes,
                    struct: rna.RNAStructure,
                    correl: pd.Series,
                    aso_coords: tuple[int, int]):
    """ Graph the base pairs in the structure as an arc plot. """
    arc_aspect = calc_arc_aspect_ratio(ax)
    for pair in struct.pairs:
        p5, p3 = pair
        center = (p5 + p3) / 2., 1.
        width = float(p3 - p5)
        height = width / arc_aspect
        arc = patches.Arc(center, width, height,
                          theta1=0., theta2=180.,
                          color=color_base_pair(pair, aso_coords, correl),
                          linewidth=0.1)
        ax.add_patch(arc)


def graph_correl(ax: plt.Axes, correl: pd.Series, aso_coords: tuple[int, int]):
    """ Graph the normalized CORR. """
    aso5, aso3 = aso_coords
    positions = pd.Series(correl.index, index=correl.index)
    where_aso = np.logical_and(positions >= aso5, positions <= aso3)
    ax.fill_between(positions, correl, 1.0,
                    where=where_aso, facecolor=ASO_COLOR)
    ax.fill_between(positions, correl, 1.0,
                    where=~where_aso, facecolor=NO_ASO_COLOR)
    ax.plot(positions[: aso5 - 1], correl[: aso5 - 1], color=NO_ASO_COLOR)
    ax.plot(positions[aso5: aso3], correl[aso5: aso3], color=ASO_COLOR)
    ax.plot(positions[aso3 + 1:], correl[aso3 + 1:], color=NO_ASO_COLOR)
    ax.plot(ax.get_xlim(), [CORR_THRESH, CORR_THRESH],
            color="#d55e00", linewidth=0.5)


def graph(struct: rna.RNAStructure,
          correl: pd.Series,
          aso_coords: tuple[int, int],
          outfile: Path):
    """ Draw the graph. """
    fig, ax = plt.subplots()
    pos_max = correl.size
    ax.set_xlim(1, pos_max)
    ax.set_ylim(-CORR_MIN, CORR_MAX + STRUCT_CORR_RATIO)
    ax.set_aspect(calc_arc_aspect_ratio(ax) / GOLDEN_RATIO)
    ax.set_xticks(np.arange(0, pos_max + 1, POSITION_STEP))
    ax.set_yticks(np.linspace(CORR_MIN, CORR_MAX, CORR_GRID))
    ax.grid(axis="x", visible=True, color="#f2f2f2", linewidth=0.25)
    ax.grid(axis="y", visible=True, color="#e5e5e5", linewidth=0.5)
    ax.spines[["left", "top", "right", "bottom"]].set_visible(False)
    plt.setp(ax.get_xticklines(), visible=False)
    plt.setp(ax.get_yticklines(), visible=False)
    ax.set_title(struct.ref)
    ax.set_xlabel("Position")
    ax.set_ylabel("Correlation")
    graph_structure(ax, struct, correl, aso_coords)
    graph_correl(ax, correl, aso_coords)
    ax.plot(aso_coords, [0., 0.], color="#f0e442")
    plt.savefig(outfile)


def main(structs_path: Path, correl_path: Path, asos_path: Path, outdir: Path):
    correl = load_correl_data(correl_path).squeeze()
    assert isinstance(correl, pd.Series)
    for i, state in enumerate(iter_rna_structs(structs_path)):
        correl_state = correl.reindex(state.table.index)
        correl_state.index = correl_state.index.get_level_values(seq.POS_NAME)
        correl_state.loc[correl_state < CORR_MIN] = CORR_MIN
        correl_state.loc[correl_state > CORR_MAX] = CORR_MAX
        coords = get_aso_coords(asos_path, get_ref_name(structs_path,
                                                        correl_path))
        outfile = outdir.joinpath(f"{state.ref}_{i}.pdf")
        graph(state, correl_state, coords, outfile)
        print(f"Graphed {outfile}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--structs", "-s",
                        help="CT file of predicted structures")
    parser.add_argument("--correl", "-d",
                        help="CSV file of CORRs")
    parser.add_argument("--asos", "-a",
                        help="CSV file of ASO coordinates for each reference")
    parser.add_argument("--outdir", "-o",
                        help="Path of the output directory")
    args = parser.parse_args()
    main(Path(args.structs),
         Path(args.correl),
         Path(args.asos),
         Path(args.outdir))
