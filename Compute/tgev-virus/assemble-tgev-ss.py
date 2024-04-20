"""
Assemble the secondary structure of the TGEV genome from the five overlapping
segments of the predicted structure.
"""

import os

import pandas as pd


def get_db_structure(db_file: str, end5: int, end3: int):
    with open(db_file) as f:
        f.readline()
        f.readline()
        dbs = f.readline().strip()
    return dict(zip(range(end5, end3 + 1), dbs, strict=True))


# Load the 5' and 3' ends of each section.
sections = pd.read_csv("sections-fold.csv", index_col="Section")
end5s = sections["5' End"]
end3s = sections["3' End"]


# Load the structure for each section.
structures = dict()
for section in sections.index:
    db_file = f"out/tgev-full-pool/fold/tgev/{section}/clip__average.db"
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
db_file = "out/tgev-full-pool/fold/tgev/full/clip__average.db"
os.mkdir(os.path.dirname(db_file))
with open(db_file, "w") as f:
    f.write(structure)

