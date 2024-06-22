from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from seismicrna.cluster.compare import assign_clusterings
from seismicrna.cluster.report import ClusterReport, NumUniqReadKeptF, NumClustsF
from seismicrna.core.header import parse_header
from seismicrna.mask.report import MaskReport, NumReadsKeptF
from seismicrna.table.base import MUTAT_REL
from seismicrna.table.load import ClustFreqTableLoader, load_pos_table

CLUSTERS = ((1, 1),
            (2, 1),
            (2, 2),
            (2, 3),
            (3, 1),
            (3, 2),
            (3, 3),
            (4, 1),
            (4, 2),
            (4, 3))
LIBRARIES = ((280, "ampl2"),
             (280, "frag2"),
             (560, "frag2"),
             (1120, "frag2"))
NUM_READS = (10000, 20000, 50000, 100000, 200000, 500000, 1000000)


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
    return Path("sim", "samples", sample, "table", ref, "full")


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


def find_report_file(step: str,
                     length: int,
                     order: int,
                     props: int,
                     library: str,
                     n_reads: int):
    """ Find the path to a report file. """
    sample = format_sample_name(order, props, library, n_reads)
    ref = format_ref_name(length)
    return Path("sim", "samples", sample, step, ref, "full", f"{step}-report.json")


def find_mask_report_file(length: int,
                          order: int,
                          props: int,
                          library: str,
                          n_reads: int):
    """ Find the path to a mask report file. """
    return find_report_file("mask", length, order, props, library, n_reads)


def find_clust_report_file(length: int,
                           order: int,
                           props: int,
                           library: str,
                           n_reads: int):
    """ Find the path to a cluster report file. """
    return find_report_file("cluster", length, order, props, library, n_reads)


def get_num_reads(mask_report: Path):
    report = MaskReport.load(mask_report)
    return int(report.get_field(NumReadsKeptF))


def get_num_clusters(cluster_report: Path):
    report = ClusterReport.load(cluster_report)
    return int(report.get_field(NumClustsF))


def get_num_uniq_reads(cluster_report: Path):
    report = ClusterReport.load(cluster_report)
    return int(report.get_field(NumUniqReadKeptF))


def iter_attrs(libraries=LIBRARIES,
               ns_reads=NUM_READS):
    """ Iterate through all combinations of attributes. """
    for length, library in libraries:
        for n_reads in ns_reads:
            yield length, library, n_reads


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


def load_observed_mus(table_file: str | Path):
    """ Load observed mutation rates from a table file. """
    table = load_pos_table(Path(table_file))
    return table.fetch_ratio(rel=MUTAT_REL)[MUTAT_REL, table.header.max_order]


def load_expected_pis(csv_file: str | Path):
    """ Load expected cluster proportions from a CSV file. """
    data = pd.read_csv(csv_file, index_col=list(range(2)))
    clusters = parse_header(data.index)
    return data.loc[clusters.max_order, "Proportion"]


def load_observed_pis(csv_file: str | Path,
                      expected_mus: pd.DataFrame,
                      observed_mus: pd.DataFrame):
    """ Load observed cluster proportions from a CSV file. """
    table = ClustFreqTableLoader(csv_file)
    data = table.data.loc[table.header.max_order]
    assignments = assign_clusterings(expected_mus.values, observed_mus.values)
    data = pd.Series(data.values[assignments], index=data.index)
    return data / data.sum()


def calc_diff_mus(expected: pd.DataFrame, observed: pd.DataFrame):
    """ Compare expected and observed mutation rates of clusters. """
    assignments = assign_clusterings(expected.values, observed.values)
    return observed.values[:, assignments] - expected.values


def calc_diff_pis(expected_pis: pd.Series,
                  observed_pis: pd.Series):
    """ Compare expected and observed mutation rates of clusters. """
    return observed_pis.values - expected_pis.values


def calc_rms(x: np.ndarray):
    """ Root-mean-square. """
    return np.sqrt(np.nanmean(np.square(x)))


def calc_norm(x: np.ndarray):
    """ L2 norm. """
    return np.sqrt(np.nansum(np.square(x)))


def calc_data_order(order: int, props: int):
    """ DataFrame of the difference between every observed and expected
    mutation rates in every simulated dataset. """
    data = list()
    for length, library, n_reads in iter_attrs():
        mask_report_file = find_mask_report_file(length,
                                                 order,
                                                 props,
                                                 library,
                                                 n_reads)
        clust_report_file = find_clust_report_file(length,
                                                   order,
                                                   props,
                                                   library,
                                                   n_reads)
        num_reads = get_num_reads(mask_report_file)
        num_uniq_reads = get_num_uniq_reads(clust_report_file)
        num_clusters = get_num_clusters(clust_report_file)
        attrs = {"ReferenceLength": length,
                 "ExpectedClusters": order,
                 "ExpectedProportions": props,
                 "NumReads": num_reads,
                 "NumUniqReads": num_uniq_reads,
                 "Library": library,
                 "ObservedClusters": num_clusters}
        if num_clusters == order:
            observed_mus = load_observed_mus(find_pos_table_file(length,
                                                                 order,
                                                                 props,
                                                                 library,
                                                                 n_reads))
            expected_mus = load_expected_mus(find_muts_param_file(length, order))
            rmsd_mus = calc_rms(calc_diff_mus(expected_mus, observed_mus))
            observed_pis = load_observed_pis(find_clust_table_file(length,
                                                                   order,
                                                                   props,
                                                                   library,
                                                                   n_reads),
                                             expected_mus,
                                             observed_mus)
            expected_pis = load_expected_pis(find_clusts_param_file(order, props))
            norm_pis = calc_norm(calc_diff_pis(expected_pis, observed_pis))
            attrs["MutationRMSD"] = rmsd_mus
            attrs["ProportionNorm"] = norm_pis
            for cluster in range(1, order + 1):
                attrs[f"Cluster {cluster}"] = observed_pis[cluster]
        data.append(attrs)
    return pd.DataFrame.from_records(data)


def graph_rmsds():
    for order, props in CLUSTERS:
        print(order, props)
        # Select data.
        data = calc_data_order(order, props)
        selector = np.logical_and.reduce([data["ExpectedClusters"] == order,
                                          data["ExpectedProportions"] == props])
        selected = data.loc[selector]
        # Number of clusters.
        fig, ax = plt.subplots()
        sns.lineplot(selected,
                     y="ObservedClusters",
                     x="NumReads",
                     hue="ReferenceLength",
                     style="Library",
                     markers=True)
        ax.set_xscale("log")
        ax.set_xlim((1.e3, 1.e6))
        ax.set_ylim((0, 5))
        plt.show()
        # Proportion bars.
        fig, ax = plt.subplots()
        bottom = None
        for cluster in range(1, order + 1):
            cname = f"Cluster {cluster}"
            x = np.arange(selected.index.size)
            y = selected.loc[:, cname]
            ax.bar(x, y, bottom=bottom, label=cname)
            bottom = bottom + y if bottom is not None else y
        xlabels = [" ".join([selected.loc[row, "Library"],
                             str(selected.loc[row, "NumReads"])])
                   for row in selected.index]
        ax.set_xlabel(xlabels)
        plt.show()
        # Proportion norm.
        fig, ax = plt.subplots()
        sns.lineplot(selected,
                     y="ProportionNorm",
                     x="NumReads",
                     hue="ReferenceLength",
                     style="Library",
                     markers=True)
        ax.set_xscale("log")
        ax.set_xlim((1.e3, 1.e6))
        plt.show()
        # Mutation RMSD.
        fig, ax = plt.subplots()
        sns.lineplot(selected,
                     y="MutationRMSD",
                     x="NumReads",
                     hue="ReferenceLength",
                     style="Library",
                     markers=True)
        ax.set_xscale("log")
        ax.set_xlim((1.e3, 1.e6))
        plt.show()


def main():
    graph_rmsds()


if __name__ == "__main__":
    main()
