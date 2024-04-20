import os
import re

from Bio import SeqIO


exp_up = 100
exp_dn = 1899


long_fse_dir = "long_fse_seqs"
slippery = re.compile("UUUAAAC")
max_deviation = 30  # nucleotides


def get_frameshift_locs(features_table):
    fse_coords = dict()
    accession = None
    fse_coord = None
    get_coord = False
    with open(features_table) as f:
        for line in f:
            if line.strip() == "":
                continue
            if line.startswith(">"):
                if accession is not None:
                    if accession not in fse_coords:
                        print("Missing", accession)
                _, accession, _ = line.split("|")
                continue
            fields = line.rstrip().split("\t")
            if len(fields) == 3 and fields[2] == "CDS":
                if fields[0].startswith("<"):
                    fse_coord = int(fields[0][1:])
                    assert accession is not None
                    assert accession not in fse_coords
                    fse_coords[accession] = fse_coord
                    accession = None
                    fse_coord = None
                else:
                    fse_coord = int(fields[1])
            elif len(fields) == 2 and fse_coord is not None:
                assert fse_coord == int(fields[0])
                assert accession is not None
                assert accession not in fse_coords
                fse_coords[accession] = fse_coord
                accession = None
                fse_coord = None
    # The FSE for TGEV (NC_038861.1) is not annotated in RefSeq. Add it manually.
    fse_coords["NC_038861.1"] = 12338
    return fse_coords


def extract_long_fse(genome_seq_file, features_file):
    if not os.path.isdir(long_fse_dir):
        os.mkdir(long_fse_dir)
    fse_coords = get_frameshift_locs(features_file)
    records = SeqIO.parse(genome_seq_file, "fasta")
    for record in records:
        accession = record.id
        seq = str(record.seq).replace("T", "U")
        fse_coord_est = fse_coords[accession]
        slippery_locs = list(slippery.finditer(seq))
        distances = [abs(fse_coord_est - x.end()) for x in slippery_locs]
        min_dist = min(distances)
        assert min_dist <= max_deviation
        closest_loc = [loc for loc, dist in zip(slippery_locs, distances) if dist == min_dist]
        assert len(closest_loc) == 1
        fse_coord = closest_loc[0].start()
        long_fse = seq[fse_coord - exp_up: fse_coord + exp_dn + 1]
        line = f">{accession}\n{long_fse}\n"
        long_fse_file = os.path.join(long_fse_dir, f"{accession}.fasta")
        with open(long_fse_file, "w") as fo:
            fo.write(line)
        
extract_long_fse("cov_refseq.fasta", "cov_features.txt")

