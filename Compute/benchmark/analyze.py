import os
import sys
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from seismicrna.cluster.compare import assign_clusterings
from seismicrna.cluster.report import ClusterReport, NumUniqReadKeptF
from seismicrna.core.header import ORDER_NAME, parse_header
from seismicrna.table.base import MUTAT_REL, UNAMB_REL
from seismicrna.table.load import load_pos_table

LENGTHS = [280, 560, 1120]
CLUSTERS = [(1, 1),
            (2, 1),
            (2, 2),
            (2, 3),
            (2, 4),
            (3, 1),
            (3, 2),
            (3, 3),
            (4, 1),
            (4, 2),
            (4, 3)]
MRATES = [1, 3, 6]
LIBRARIES = ["frag1", "frag2", "ampl2"]
NUM_READS = [10000]


def format_sample_name(order: int,
                       props: int,
                       mrate: int,
                       library: str,
                       n_reads: int):
    return f"c{order}-{props}-m{mrate}-{library}-n{n_reads}"


def format_ref_name(length: int):
    return f"ref-{length}"


def find_muts_param_file(length: int, order: int, mrate: int):
    """ Find the path to a mutation rate parameter file. """
    filename = f"c{order}-m{mrate}.muts.csv"
    ref = format_ref_name(length)
    return Path("sim", "params", ref, "full", filename)


def find_clusts_param_file(order: int, proportions: int):
    """ Find the path to a cluster proportions parameter file. """
    filename = f"c{order}-{proportions}.csv"
    return Path("clusts", filename)


def find_table_dir(length: int,
                   order: int,
                   props: int,
                   mrate: int,
                   library: str,
                   n_reads: int):
    """ Find the directory of a table file. """
    sample = format_sample_name(order, props, mrate, library, n_reads)
    ref = format_ref_name(length)
    return Path("out", sample, "table", ref, "full")


def find_pos_table_file(length: int,
                        order: int,
                        props: int,
                        mrate: int,
                        library: str,
                        n_reads: int):
    """ Find the path to a table file of positions. """
    return Path(find_table_dir(length, order, props, mrate, library, n_reads),
                "clust-per-pos.csv")


def find_clust_table_file(length: int,
                          order: int,
                          props: int,
                          mrate: int,
                          library: str,
                          n_reads: int):
    """ Find the path to a table file of cluster proportions. """
    return Path(find_table_dir(length, order, props, mrate, library, n_reads),
                "clust-freq.csv")


def find_clust_report_file(length: int,
                           order: int,
                           props: int,
                           mrate: int,
                           library: str,
                           n_reads: int):
    """ Find the path to a table file of cluster proportions. """
    sample = format_sample_name(order, props, mrate, library, n_reads)
    ref = format_ref_name(length)
    return Path("out", sample, "cluster", ref, "full", "cluster-report.json")


def get_num_uniq_reads(cluster_report: Path):
    report = ClusterReport.load(cluster_report)
    return int(report.get_field(NumUniqReadKeptF))


def iter_attrs():
    """ Iterate through all combinations of attributes. """
    for length in LENGTHS:
        for order, props in CLUSTERS:
            for mrate in MRATES:
                for library in LIBRARIES:
                    for n_reads in NUM_READS:
                        yield length, order, props, mrate, library, n_reads


def load_expected(csv_file: str | Path):
    """ Load expected mutation rates from a CSV file. """
    data = pd.read_csv(csv_file,
                       index_col=list(range(2)),
                       header=list(range(3)))
    # Cast the columns from str to int.
    data.columns = parse_header(data.columns).index
    raw_mut_rate = (data.loc[:, "16"]
                    + data.loc[:, "32"]
                    + data.loc[:, "64"]
                    + data.loc[:, "128"])
    return raw_mut_rate / (raw_mut_rate + data.loc[:, "1"])


def load_observed(table_file: str | Path):
    """ Load observed mutation rates from a table file. """
    table = load_pos_table(Path(table_file))
    return table.fetch_ratio(rel=MUTAT_REL)[MUTAT_REL]


def load_coverage(table_file: str | Path):
    """ Load read coverage from a table file. """
    table = load_pos_table(Path(table_file))
    return table.fetch_count(rel=UNAMB_REL)[UNAMB_REL]


def calc_diffs_clusters(expected: pd.DataFrame,
                        observed: pd.DataFrame,
                        coverage: pd.DataFrame):
    """ Compare expected and observed mutation rates of clusters. """
    orders = list(set(expected.columns.get_level_values(ORDER_NAME).to_list()))
    if len(orders) != 1:
        raise ValueError(f"Expected one order, but got\n{expected}")
    order = int(orders[0])
    expected_order = expected[order].values
    try:
        observed_order = observed[order].values
        coverage_order = coverage[order].values
    except KeyError:
        return np.array([])
    assignments = assign_clusterings(expected_order, observed_order)
    expected_total = np.concatenate([expected_order[:, i]
                                     for i, _ in assignments],
                                    axis=None)
    observed_total = np.concatenate([observed_order[:, j]
                                     for _, j in assignments],
                                    axis=None)
    coverage_total = np.concatenate([coverage_order[:, j]
                                     for _, j in assignments],
                                    axis=None)
    diffs = observed_total - expected_total
    mask = ~np.isnan(diffs)
    return coverage_total[mask], diffs[mask]


def calc_diffs():
    """ DataFrame of the difference between every observed and expected
    mutation rates in every simulated dataset. """
    data = list()
    for length, order, props, mrate, library, n_reads in iter_attrs():
        table_file = find_pos_table_file(length,
                                         order,
                                         props,
                                         mrate,
                                         library,
                                         n_reads)
        if not table_file.is_file():
            continue
        observed = load_observed(table_file)
        coverage = load_coverage(table_file)
        expected = load_expected(find_muts_param_file(length, order, mrate))
        for cover, diff in zip(*calc_diffs_clusters(expected,
                                                    observed,
                                                    coverage)):
            data.append({"ReferenceLength": length,
                         "Clusters": order,
                         "Proportions": props,
                         "MutationRate": mrate,
                         "Library": library,
                         "Coverage": cover,
                         "Difference": diff})
    return pd.DataFrame.from_records(data)


def graph_diffs():
    diffs = calc_diffs()
    for length in LENGTHS:
        for order, props in CLUSTERS:
            data = diffs.loc[np.logical_and(diffs["ReferenceLength"] == length,
                                            diffs["Clusters"] == order,
                                            diffs["Proportions"] == props)]
            print(length, order, props)
            sns.scatterplot(data, y="Difference", x="Coverage", hue="Library")
            plt.show()


def main():
    graph_diffs()


if __name__ == "__main__":
    main()
