from collections import defaultdict
from pathlib import Path
from typing import Iterable

from seismicrna.core.rna import from_ct


PAIRS = {
    "LS1": {
        (895, 1867),
        (896, 1866),
        (897, 1865),
        (898, 1864),
        (899, 1863),
        (900, 1862),
        (901, 1861),
        (903, 1859),
        (904, 1858),
        (905, 1857),
        (906, 1856),
        (907, 1855),
        (908, 1854),
        (909, 1853),
    },
    "LS2a": {
        (872, 1978),
        (873, 1977),
        (874, 1976),
        (876, 1975),
        (877, 1974),
        (878, 1973),
        (879, 1972),
        (880, 1971),
        (881, 1970),
        (882, 1969),
        (883, 1968),
        (886, 1964),
        (887, 1963),
        (888, 1962),
    },
    "LS2b": {
        (849, 1995),
        (850, 1994),
        (851, 1993),
        (853, 1991),
        (854, 1990),
        (855, 1989),
        (858, 1988),
        (859, 1987),
        (860, 1986),
        (861, 1985),
        (862, 1984),
        (863, 1983),
        (864, 1982),
    },
}
PAIRS["all"] = PAIRS["LS1"] | PAIRS["LS2a"] | PAIRS["LS2b"]

FOLD_DIR = Path("out/ctrl2-deep/fold/sars2-2924/sars2-1799")
CT_FILES = [
    "fse9__cluster-1-1.ct",
    "fse9__cluster-2-1.ct",
    "fse9__cluster-2-2.ct",
    "fse9__nodms.ct",
]


def check_pairs(query: Iterable[tuple[int, int]],
                truth: Iterable[tuple[int, int]]):
    return not bool(set(truth) - set(query))


def check_file(ct_file: Path):
    print(f"Checking {ct_file}")
    long_range = defaultdict(list)
    for structure in from_ct(ct_file):
        for model, pairs in PAIRS.items():
            long_range[model].append(check_pairs(structure.pairs, pairs))
    for model, structures in long_range.items():
        fraction = sum(structures) / len(structures)
        print(f"{model}:\t{fraction}")


def check_files():
    for ct_file in CT_FILES:
        check_file(FOLD_DIR.joinpath(ct_file))


if __name__ == "__main__":
    check_files()

