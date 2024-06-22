# Benchmark the accurary of SEISMIC-RNA

import os
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

mean_fragment_length = 200
nprofiles = {1: 1, 2: 3, 3: 3, 4: 3}


def link_if_needed(src: str, dst: str):
    try:
        os.symlink(src, dst)
    except FileExistsError:
        pass


def simulate():
    """ Simulate samples for each reference sequence. """
    relate_reports = list()
    for length in [280, 560, 1120]:
        ref = f"ref-{length}"
        pdir = os.path.join(os.getcwd(), "sim", "params", ref, "full")
        param3 = {"frag2": mean_fragment_length / (length * 2.) + 0.5, "ampl2": 1.0}
        paraml = {"frag2": mean_fragment_length / length, "ampl2": 1.0}
        fasta = f"sim/refs/{ref}.fa"
        print("ref")
        ref_mod.run(refs=ref, ref=ref, reflen=length)
        for clusters in range(1, MAX_CLUSTERS + 1):
            cname = f"c{clusters}"
            ct = os.path.join(pdir, f"{cname}.ct")
            print("fold")
            fold_mod.run(fasta=fasta, profile_name=cname, fold_max=clusters)
            muts = os.path.join(pdir, f"{cname}.muts.csv")
            print("muts")
            muts_mod.run(ct_file=(ct,))
            for library in ["ampl2", "frag2"]:
                if library == "ampl2" and length != 280:
                    continue
                lct = os.path.join(pdir, f"{library}.ct")
                link_if_needed(ct, lct)
                ends = os.path.join(pdir, f"{library}.ends.csv")
                ends_mod.run(ct_file=(lct,),
                             end3_fmean=param3[library],
                             insert_fmean=paraml[library],
                             ends_var=ENDS_VAR)
                for profile in range(1, nprofiles[clusters] + 1):
                    pname = f"c{clusters}-{profile}-{library}"
                    pct = os.path.join(pdir, f"{pname}.ct")
                    link_if_needed(ct, pct)
                    pends = os.path.join(pdir, f"{pname}.ends.csv")
                    link_if_needed(ends, pends)
                    pmuts = os.path.join(pdir, f"{pname}.muts.csv")
                    link_if_needed(muts, pmuts)
                    pclusts = os.path.join(pdir, f"{pname}.clusts.csv")
                    link_if_needed(os.path.join(os.getcwd(), "clusts", f"c{clusters}-{profile}.csv"),
                                   pclusts)
                    for reads in [10000, 20000, 50000, 100000, 200000, 500000, 1000000]:
                        sample = f"{pname}-n{reads}"
                        print(f"Sample {sample} for ref {ref}")
                        relate_reports.extend(relate_mod.run(param_dir=(pdir,),
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
