import os

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from seismicrna.core.rna import from_ct

# Editable text in PDF.
plt.rcParams["pdf.fonttype"] = 42

fse_structs_dir = "long_fse_structs"
seq_len = 2000


def get_ct_files():
    return [os.path.join(fse_structs_dir, f)
            for f in os.listdir(fse_structs_dir)]


def count_pairs(structures, range_start=101, range_end=250):
    n = len(structures[0].seq)
    counts = np.zeros(n, dtype=int)
    for structure in structures:
        for base5, base3 in structure.pairs:
            assert 0 < base5 < base3 <= n
            in5 = range_start <= base5 <= range_end
            in3 = range_start <= base3 <= range_end
            if in5 and not in3:
                counts[base3 - 1] += 1
            elif in3 and not in5:
                counts[base5 - 1] += 1
    return counts


def get_contact_freqs():
    ct_files = get_ct_files()
    accs = [f[f.index("/") + 1: -3] for f in ct_files]
    pos = list(range(1, seq_len + 1))
    contact_freqs = pd.DataFrame(
        np.zeros((len(accs), seq_len), dtype=float),
        index=accs,
        columns=pos,
    )
    n_structures = dict()
    for acc, ct_file in zip(accs, ct_files):
        print(acc)
        structures = list(from_ct(ct_file))
        n_structures[acc] = len(structures)
        contact_freqs.loc[acc] = count_pairs(structures) / n_structures[acc]
    return contact_freqs, n_structures


def plot_contact_freqs(contact_freqs):
    contact_freqs.to_csv("contact_freqs.csv")
    sns.clustermap(contact_freqs, row_cluster=True, col_cluster=False, yticklabels=contact_freqs.index)
    plt.savefig("contact_freqs.pdf")
    plt.close()

f, n = get_contact_freqs()
print("min structures", min(n.values()))
print("max structures", max(n.values()))
plot_contact_freqs(f)

