from typing import Iterable

import numpy as np
import pandas as pd


LRI_FILE = "models/sars2/full/fse9__cluster-2-1.db"
PK_FILE = "models/sars2/pk/zhang2021.db"
MIX_FILE = "models/sars2/mixed/mixed.db"
PARTS_FILE = "model-parts.csv"


COORDS = {
    "AS1": ((252, 261), (298, 307)),
    "LS2b": ((363, 369), (1503, 1509)),
    "LS3": ((329, 359), (1521, 1550)),
    "LS3-trunc": ((336, 359), (1521, 1543)),
    "LS4": ((309, 326), (1600, 1622)),
    "PS1": ((305, 316), (322, 332)),
    "PS2": ((317, 321), (367, 371)),
    "PS3": ((334, 342), (352, 361)),
}


def read_db_struct(file: str, pad5: int = 0, pad3: int = 0):
    with open(file) as f:
        f.readline()
        seq = f.readline().rstrip()
        return seq, np.hstack([
            np.full(pad5, "."),
            np.array(list(f.readline().rstrip())),
            np.full(pad3, "."),
        ])


def combine_model(struct1: np.ndarray,
                  parts1: Iterable[str],
                  struct2: np.ndarray,
                  parts2: Iterable[str]):
    if (size := struct1.size) != struct2.size:
        raise ValueError(f"Sizes disagree between struct1 ({struct1.size}) "
                         f"and struct2 ({struct2.size}")
    combined = np.full(size, ".")
    for struct, parts in [(struct1, parts1), (struct2, parts2)]:
        for part in parts:
            for start, end in COORDS[part]:
                combined[start - 1: end] = struct[start - 1: end]
    return combined


def make_model(lri: np.ndarray, pk: np.ndarray, parts: list[str]):
    lri_parts = [part for part in parts if not part.startswith("P")]
    pk_parts = [part for part in parts if part.startswith("P")]
    return combine_model(lri, lri_parts, pk, pk_parts)


def write_db_structs(file: str, seq: str, models: dict[str, np.ndarray]):
    with open(file, "x") as f:
        for i, (model, struct) in enumerate(models.items()):
            f.write(f">{model}\n")
            if i == 0:
                f.write(f"{seq}\n")
            f.write(f"{''.join(struct)}\n")


def main():
    seq, lri = read_db_struct(LRI_FILE)
    _, pk = read_db_struct(PK_FILE, pad5=302, pad3=1428)
    mixed = dict()
    for model, parts in pd.read_csv(PARTS_FILE, index_col="Model").T.items():
        mixed[model] = make_model(lri, pk, parts.loc[parts == 1].index.to_list())
    write_db_structs(MIX_FILE, seq, mixed)


if __name__ == "__main__":
    main()

