""" Copy references in use from all_refs.bib to refs.bib """

import re
from string import ascii_letters, digits


NAME_CHARS = set(ascii_letters + digits)
BBL_NAME_KEY = r"\bibitem{"

BBL_FILE = "main.bbl"
FROM_BIB_FILE = "all_refs.bib"
TO_BIB_FILE = "refs.bib"


def parse_bib_name(line: str):
    bracket = line.index("{")
    comma = line.index(",")
    name = line[bracket + 1: comma]
    if extra := (set(name) - NAME_CHARS):
        raise ValueError(f"Name {repr(name)} has invalid characters: {extra}")
    return name


def parse_bib(bib_file: str):
    with open(bib_file) as f:
        entry = list()
        for i, line in enumerate(f, start=1):
            rs = line.rstrip()
            if not rs:
                # Skip blank lines.
                continue
            if rs.startswith("@"):
                # The line begins an entry.
                if entry:
                    raise ValueError(f"Line {i}: entry {rs} began before entry "
                                     f"{entry} ended")
            elif not entry:
                raise ValueError(f"Line {i}: {rs} not part of any entry")
            entry.append(rs)
            if rs == "}":
                # The line ends an entry.
                yield parse_bib_name(entry[0]), "\n".join(entry)
                entry.clear()


def parse_bbl_name(line: str):
    if not line.startswith(BBL_NAME_KEY):
        raise ValueError(
                f"Line {repr(line)} does not start with {repr(BBL_NAME_KEY)}"
        )
    bracket = line.index("}")
    name = line[len(BBL_NAME_KEY): bracket]
    if extra := (set(name) - NAME_CHARS):
        raise ValueError(f"Name {repr(name)} has invalid characters: {extra}")
    return name


def parse_bbl(bbl_file: str):
    with open(bbl_file) as f:
        for line in f:
            if line.startswith(BBL_NAME_KEY):
                yield parse_bbl_name(line)


if __name__ == "__main__":
    entries = dict(parse_bib(FROM_BIB_FILE))
    with open(TO_BIB_FILE, "x") as f:
        for name in parse_bbl(BBL_FILE):
            f.write(f"{entries[name]}\n")

