"""
Assemble the secondary structure of the TGEV genome from the five overlapping
segments of the predicted structure.
"""

import os

import pandas as pd

from seismicrna.core.seq import DNA, parse_fasta

PROFILE = "full__average"


def get_db_structure(db_file: str, end5: int, end3: int):
    with open(db_file) as f:
        f.readline()
        f.readline()
        dbs = f.readline().strip()
    return dict(zip(range(end5, end3 + 1), dbs, strict=True))


# Load the reference sequence.
refseq = dict(parse_fasta("tgev_consensus.fa", DNA))["tgev"].tr()


# Load the 5' and 3' ends of each section.
sections = pd.read_csv("sections-fold.csv", index_col="Section")
end5s = sections["5' End"]
end3s = sections["3' End"]


# Load the structure for each section.
structures = dict()
for section in sections.index:
    db_file = f"out/tgev-full-pool/fold/tgev/{section}/{PROFILE}.db"
    structures[section] = get_db_structure(db_file,
                                           end5s[section],
                                           end3s[section])


# Determine which structure to use for each position.
cutoffs = (end5s.values[1:] + end3s.values[:-1]) / 2.
def get_section_for_pos(pos: int):
    for section, cutoff in zip(sections.index[:-1], cutoffs, strict=True):
        if pos < cutoff:
            return section
    return sections.index[-1]
num_pos = end3s.max()
pos_sections = {pos: get_section_for_pos(pos) for pos in range(1, num_pos + 1)}


# Assemble the consensus structure.
structure = "".join(structures[section][pos]
                    for pos, section in pos_sections.items())

# Write the full-genome structure.
db_file = "out/tgev-full-pool/fold/tgev/full/full__average.db"
try:
    os.mkdir(os.path.dirname(db_file))
except FileExistsError:
    pass
with open(db_file, "w") as f:
    f.write(f">{PROFILE}\n")
    f.write(f"{refseq}\n")
    f.write(structure)
ct_file = "out/tgev-full-pool/fold/tgev/full/full__average.ct"
dot2ct_cmd = f"dot2ct {db_file} {ct_file}"
if status := os.system(dot2ct_cmd):
    raise ValueError(f"Command {dot2ct_cmd} exited with status {status}")

# Write the 5' UTR structure.
utr5_start = 1
utr5_end = 330
db_file = "out/tgev-full-pool/fold/tgev/utr5/full__average.db"
try:
    os.mkdir(os.path.dirname(db_file))
except FileExistsError:
    pass
with open(db_file, "w") as f:
    f.write(f">{PROFILE}\n")
    f.write(f"{refseq[utr5_start-1:utr5_end]}\n")
    f.write(structure[utr5_start-1:utr5_end])
ct_file = "out/tgev-full-pool/fold/tgev/utr5/full__average.ct"
dot2ct_cmd = f"dot2ct {db_file} {ct_file}"
if status := os.system(dot2ct_cmd):
    raise ValueError(f"Command {dot2ct_cmd} exited with status {status}")

