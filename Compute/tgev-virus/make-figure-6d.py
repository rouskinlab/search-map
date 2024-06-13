import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from matplotlib.patches import Arc
from seismicrna.core.rna import from_ct


AUCROC_FILE = "out/tgev-full-pool/graph/tgev/full/aucroll_full__masked_m-ratio-q0_45-9.csv"
STRUCT_FILE = "out/tgev-full-pool/fold/tgev/full/full__average.ct"
GRAPH_FILE = "../../MainFigures/tgev/figure-6d.svg"
END5 = 12042
END3 = 13840
MAX_WIDTH = 300
GOLDEN_RATIO = (1. + 5.**0.5) / 2.
STRUCTURE_RATIO = 1. / 3.
STRUCTURE_HEIGHT = STRUCTURE_RATIO / MAX_WIDTH
XTICK_INTERVAL = 500


def graph_aucroc(ax):
    aucroc = pd.read_csv(AUCROC_FILE,
                         header=[0, 1],
                         index_col=[0, 1]).squeeze()
    assert isinstance(aucroc, pd.Series)
    aucroc = aucroc.loc[END5: END3]
    positions = aucroc.index.get_level_values("Position")
    ax.fill_between(positions, aucroc, 1.)
    ax.plot(positions, aucroc)


def graph_struct(ax):
    structs = list(from_ct(STRUCT_FILE))
    if len(structs) != 1:
        raise ValueError(f"Expected 1 structure, but got {structs}")
    struct = structs[0]
    for end5, end3 in struct.pairs:
        if END5 <= end5 <= end3 <= END3:
            center = ((end5 + end3) / 2., 1.)
            width = end3 - end5
            if width > MAX_WIDTH:
                raise ValueError(f"Width {width} exceeds maximum {MAX_WIDTH}")
            height = STRUCTURE_HEIGHT * width
            ax.add_patch(Arc(
                center, width, height,
                theta1=0., theta2=180.,
                color="#e5e5e5",
            ))


def main():
    fig, ax = plt.subplots()
    graph_aucroc(ax)
    graph_struct(ax)
    ax.set_xlim((END5, END3))
    ax.set_ylim((0., 1. + STRUCTURE_RATIO))
    ax.set_aspect(1. / (GOLDEN_RATIO * STRUCTURE_HEIGHT))
    ax.set_xticks(np.arange(XTICK_INTERVAL * (1 + int(END5 / XTICK_INTERVAL)),
                            END3,
                            XTICK_INTERVAL))
    ax.set_yticks(np.linspace(0., 1., 5))
    plt.tight_layout()
    plt.savefig(GRAPH_FILE)
    plt.close()


if __name__ == "__main__":
    main()

