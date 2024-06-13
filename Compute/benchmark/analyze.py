import os
import sys
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from seismicrna.cluster.compare import assign_clusterings
from seismicrna.cluster.report import ClusterReport, NumUniqReadKeptF, NumClustsF
from seismicrna.core.header import ORDER_NAME, parse_header
from seismicrna.table.base import MUTAT_REL, UNAMB_REL
from seismicrna.table.load import load_pos_table

LENGTHS = (280, 560, 1120)
CLUSTERS = ((1, 1),
            (2, 1),
            (2, 2),
            (2, 3),
            (2, 4),
            (3, 1),
            (3, 2),
            (3, 3),
            (4, 1),
            (4, 2),
            (4, 3))
LIBRARIES = ("frag2", "ampl2")
NUM_READS = (10000, 30000, 100000, 300000, 1000000)


def format_sample_name(order: int,
                       props: int,
                       library: str,
                       n_reads: int):
    return f"c{order}-{props}-{library}-n{n_reads}"


def format_ref_name(length: int):
    return f"ref-{length}"


def find_muts_param_file(length: int, order: int):
    """ Find the path to a mutation rate parameter file. """
    filename = f"c{order}.muts.csv"
    ref = format_ref_name(length)
    return Path("sim", "params", ref, "full", filename)


def find_clusts_param_file(order: int, proportions: int):
    """ Find the path to a cluster proportions parameter file. """
    filename = f"c{order}-{proportions}.csv"
    return Path("clusts", filename)


def find_table_dir(length: int,
                   order: int,
                   props: int,
                   library: str,
                   n_reads: int):
    """ Find the directory of a table file. """
    sample = format_sample_name(order, props, library, n_reads)
    ref = format_ref_name(length)
    return Path("out", sample, "table", ref, "full")


def find_pos_table_file(length: int,
                        order: int,
                        props: int,
                        library: str,
                        n_reads: int):
    """ Find the path to a table file of positions. """
    return Path(find_table_dir(length, order, props, library, n_reads),
                "clust-per-pos.csv")


def find_clust_table_file(length: int,
                          order: int,
                          props: int,
                          library: str,
                          n_reads: int):
    """ Find the path to a table file of cluster proportions. """
    return Path(find_table_dir(length, order, props, library, n_reads),
                "clust-freq.csv")


def find_clust_report_file(length: int,
                           order: int,
                           props: int,
                           library: str,
                           n_reads: int):
    """ Find the path to a table file of cluster proportions. """
    sample = format_sample_name(order, props, library, n_reads)
    ref = format_ref_name(length)
    return Path("out", sample, "cluster", ref, "full", "cluster-report.json")


def get_num_clusters(cluster_report: Path):
    report = ClusterReport.load(cluster_report)
    return int(report.get_field(NumClustsF))


def get_num_uniq_reads(cluster_report: Path):
    report = ClusterReport.load(cluster_report)
    return int(report.get_field(NumUniqReadKeptF))


def iter_attrs(lengths=LENGTHS,
               clusters=CLUSTERS,
               libraries=LIBRARIES,
               ns_reads=NUM_READS):
    """ Iterate through all combinations of attributes. """
    for length in lengths:
        for order, props in clusters:
            for library in libraries:
                for n_reads in ns_reads:
                    yield length, order, props, library, n_reads


def load_expected_mus(csv_file: str | Path):
    """ Load expected mutation rates from a CSV file. """
    data = pd.read_csv(csv_file,
                       index_col=list(range(2)),
                       header=list(range(3)))
    # Cast the columns from str to int.
    clusters = parse_header(data.columns)
    data.columns = clusters.index
    raw_mut_rate = (data.loc[:, "16"]
                    + data.loc[:, "32"]
                    + data.loc[:, "64"]
                    + data.loc[:, "128"])
    mus = raw_mut_rate / (raw_mut_rate + data.loc[:, "1"])
    return mus.loc[:, clusters.max_order]


def load_observed_mus(table_file: str | Path, order: int | None = None):
    """ Load observed mutation rates from a table file. """
    table = load_pos_table(Path(table_file))
    if order is None:
        order = parse_header(table.header).max_order
    return table.fetch_ratio(rel=MUTAT_REL)[MUTAT_REL, order]


def load_expected_pis(csv_file: str | Path):
    """ Load expected cluster proportions from a CSV file. """
    data = pd.read_csv(csv_file, index_col=list(range(2)))
    clusters = parse_header(data.index)
    return data.loc[clusters.max_order, "Proportion"]


def load_observed_pis(csv_file: str | Path):
    """ Load observed cluster proportions from a CSV file. """
    data = pd.read_csv(csv_file, index_col=list(range(2)))
    clusters = parse_header(data.index)
    counts = data.loc[clusters.max_order, "Number of Reads"]
    return counts / counts.sum()


def calc_diff_mus(expected: pd.DataFrame, observed: pd.DataFrame):
    """ Compare expected and observed mutation rates of clusters. """
    assignments = assign_clusterings(expected.values, observed.values)
    return observed.values[:, assignments] - expected.values


def calc_rms(x: np.ndarray):
    """ Root-mean-square. """
    return np.sqrt(np.nanmean(np.square(x)))


def calc_data():
    """ DataFrame of the difference between every observed and expected
    mutation rates in every simulated dataset. """
    nums_clusters = list()
    rmsds_mus = list()
    rmsds_pis = list()
    for length, order, props, library, n_reads in iter_attrs(lengths=[280]):
        clust_report_file = find_clust_report_file(length,
                                                   order,
                                                   props,
                                                   library,
                                                   n_reads)
        num_uniq_reads = get_num_uniq_reads(clust_report_file)
        num_clusters = get_num_clusters(clust_report_file)
        attrs = {"ReferenceLength": length,
                 "ExpectedClusters": order,
                 "ExpectedProportions": props,
                 "NumUniqReads": num_uniq_reads,
                 "Library": library}
        nums_clusters.append(attrs | {"ObservedClusters": num_clusters})
        if num_clusters != order:
            continue
        observed_mus = load_observed_mus(find_pos_table_file(length,
                                                             order,
                                                             props,
                                                             library,
                                                             n_reads),
                                         order)
        expected_mus = load_expected_mus(find_muts_param_file(length, order))
        rmsd_mus = calc_rms(calc_diff_mus(expected_mus, observed_mus))
        rmsds_mus.append(attrs | {"MutationRMSD": rmsd_mus})
    nums_clusters = pd.DataFrame.from_records(nums_clusters)
    rmsds_mus = pd.DataFrame.from_records(rmsds_mus)
    rmsds_pis = pd.DataFrame.from_records(rmsds_pis)
    return nums_clusters, rmsds_mus, rmsds_pis


def graph_rmsds():
    nums_clusters, rmsds_mus, rmsds_pis = calc_data()
    for length in [280]:
        for clusts, props in CLUSTERS:
            print(clusts, props)
            select = np.logical_and.reduce([nums_clusters["ReferenceLength"] == length,
                                            nums_clusters["ExpectedClusters"] == clusts,
                                            nums_clusters["ExpectedProportions"] == props])
            data = nums_clusters.loc[select]
            fig, ax = plt.subplots()
            sns.lineplot(data,
                         y="ObservedClusters",
                         x="NumUniqReads",
                         hue="ReferenceLength",
                         style="Library")
            ax.set_ylim((0, 6))
            plt.show()


def main():
    graph_rmsds()


if __name__ == "__main__":
    main()
