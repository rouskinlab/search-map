#!python

import os
import sys

from seismicrna.relate.aux.cigarop import count_cigar_ref

SEP = "\t"
HEADER = "@"
MIN_FIELDS = 11
NAME_FIELD = 0
CIGAR_FIELD = 5
POS_FIELD = 3
REV_TOLERANCE = 2
FWD_TOLERANCE = 6
POSITIONS = [
    (657, 939),  # FSE Amplicon
    (1793, 2082),  # Amplicon 9
]


def compute_span(line: str):
    """ Get the coordinates of the reference spanned by a SAM line. """
    fields = line.split(SEP)
    if len(fields) < MIN_FIELDS:
        raise ValueError(f"Invalid SAM line: {repr(line)}")
    pos5 = int(fields[POS_FIELD])
    cigar = fields[CIGAR_FIELD]
    pos3 = pos5 + count_cigar_ref(cigar) - 1
    return pos5, pos3


def pair_passes(line1: str, line2: str, valid5: set[int], valid3: set[int]):
    """ Whether the pair maps to a valid span of positions. """
    span15, span13 = compute_span(line1)
    span25, span23 = compute_span(line2)
    if span15 in valid5 and span23 in valid3:
        # Line 1 is upstream.
        return True
    if span13 in valid3 and span25 in valid5:
        # Line 2 is upstream.
        return True
    return False


def filter_sam(input_file: str,
               output_file: str,
               positions: list[tuple[int, int]],
               out_tol: int = 0,
               in_tol: int = 0):
    """ Filter a SAM file, keeping reads mapping no more than `rev_tol`
    behind and `fwd_tol` ahead of any position in `positions`. """
    assert out_tol >= 0
    assert in_tol >= 0
    valid5 = set()
    valid3 = set()
    for offset in range(-out_tol, in_tol + 1):
        for pos5, pos3 in positions:
            valid5.add(pos5 + offset)
            valid3.add(pos3 - offset)
    with open(input_file) as fi, open(output_file, "w") as fo:
        n_written = 0
        n_skipped = 0
        # Read every line.
        for line in fi:
            if not line.startswith(HEADER):
                # Also read the line of this line's mate.
                line2 = fi.readline()
                if ((name1 := line.split(SEP)[NAME_FIELD])
                        != (name2 := line2.split(SEP)[NAME_FIELD])):
                    raise ValueError(f"Names differ: {name1} ≠ {name2}")
                if pair_passes(line, line2, valid5, valid3):
                    # The line passed the filter, so write it to output.
                    fo.write(line)
                    # Also write the mate line.
                    fo.write(line2)
                    n_written += 1
                else:
                    n_skipped += 1
            else:
                # Write every header line.
                fo.write(line)
        try:
            f_written = n_written / (n_written + n_skipped)
        except ZeroDivisionError:
            f_written = 0.
        print(f"Wrote {n_written} pair(s) ({round(f_written * 100, 1)} %, "
              f"skipped {n_skipped} pairs) mapping to {sorted(valid5)}, "
              f"{sorted(valid3)} from {input_file} to {output_file}")


def main(input_file: str, output_file: str):
    filter_sam(input_file, output_file, POSITIONS, REV_TOLERANCE, FWD_TOLERANCE)


if __name__ == "__main__":
    main(*sys.argv[1:])
