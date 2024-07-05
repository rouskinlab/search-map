# Measure the accurary of SEISMIC-RNA

import os
from itertools import chain
from pathlib import Path

from seismicrna.sim import (ends as ends_mod,
                            fold as fold_mod,
                            muts as muts_mod,
                            relate as relate_mod,
                            ref as ref_mod)
from seismicrna import (mask as mask_mod,
                        cluster as cluster_mod,
                        table as table_mod)
from seismicrna.core import logs

logs.set_config(2, 0)

MAX_CLUSTERS = 4
ENDS_VAR = 1. / 3.
LENGTHS = [280, 560, 1120]
SMALLEST = LENGTHS[0]
REPS = 12

mean_fragment_length = 200
nprofiles = {1: 1, 2: 3, 3: 3, 4: 3}


def link_if_needed(src: str, dst: str):
    try:
        os.symlink(os.path.abspath(src), dst)
    except FileExistsError:
        pass


def simulate():
    """ Simulate samples for each reference sequence. """
    relate_reports = list()
    for length in LENGTHS:
        refs = [f"ref-{length}"]
        if length == SMALLEST:
            refs.extend(f"ref-{length}-{i}" for i in range(1, REPS))
        pdirs = [os.path.join(os.getcwd(), "sim", "params", ref, "full")
                 for ref in refs]
        param3 = {"frag2": mean_fragment_length / (length * 2.) + 0.5, "ampl2": 1.0}
        paraml = {"frag2": mean_fragment_length / length, "ampl2": 1.0}
        print("ref")
        fastas = [ref_mod.run(refs=ref, ref=ref, reflen=length) for ref in refs]
        for clusters in range(1, MAX_CLUSTERS + 1):
            cname = f"c{clusters}"
            print("fold")
            cts = list(chain(*[fold_mod.run(fasta=fasta, profile_name=cname, fold_max=clusters)
                               for fasta in fastas]))
            pdirs = [os.path.dirname(ct) for ct in cts]
            print("PDIRS")
            print(pdirs)
            print("muts")
            muts = muts_mod.run(ct_file=cts)
            for library in ["ampl2", "frag2"]:
                if library == "ampl2" and length != SMALLEST:
                    continue
                lcts = list()
                for pdir, ct in zip(pdirs, cts, strict=True):
                    lct = os.path.join(pdir, f"{library}.ct")
                    lcts.append(lct)
                    link_if_needed(ct, lct)
                print("LCTS")
                print(lcts)
                ends = ends_mod.run(ct_file=lcts,
                                    end3_fmean=param3[library],
                                    insert_fmean=paraml[library],
                                    ends_var=ENDS_VAR)
                for profile in range(1, nprofiles[clusters] + 1):
                    pname = f"c{clusters}-{profile}-{library}"
                    print(len(pdirs), len(cts), len(ends), len(muts))
                    for pdir, ct, e, m in zip(pdirs, cts, ends, muts, strict=True):
                        pct = os.path.join(pdir, f"{pname}.ct")
                        link_if_needed(ct, pct)
                        pends = os.path.join(pdir, f"{pname}.ends.csv")
                        link_if_needed(e, pends)
                        pmuts = os.path.join(pdir, f"{pname}.muts.csv")
                        link_if_needed(m, pmuts)
                        pclusts = os.path.join(pdir, f"{pname}.clusts.csv")
                        link_if_needed(os.path.join(os.getcwd(), "clusts", f"c{clusters}-{profile}.csv"),
                                       pclusts)
                    for reads in [2000, 3000, 4000,
                                  5000, 10000, 20000,
                                  50000, 100000, 200000,
                                  500000, 1000000]:
                        sample = f"{pname}-n{reads}"
                        print(f"Sample {sample} for ref {refs}")
                        relate_reports.extend(relate_mod.run(param_dir=pdirs,
                                                             profile_name=pname,
                                                             sample=sample,
                                                             num_reads=reads,
                                                             paired_end=True))
    return relate_reports


def process(relate_reports: list[Path]):
    """ Process the simulated samples with SEISMIC-RNA. """
    mask_reports = mask_mod.run(input_path=relate_reports,
                                mask_sections_file=None,
                                mask_pos_file=None)
    cluster_reports = cluster_mod.run(input_path=mask_reports,
                                      max_clusters=(MAX_CLUSTERS + 1))
    tables = table_mod.run(input_path=cluster_reports, table_read=False)
    return tables


if __name__ == "__main__":
    if not process(simulate()):
        raise ValueError("No tables were generated")
