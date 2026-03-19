# Copyright (c) 2019 10x Genomics, Inc. All rights reserved.
"""
Functions for performing joint cell calling.
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import pandas as pd
from sklearn.neighbors import NearestCentroid
from sklearn.cluster import KMeans
from six import ensure_binary, iteritems, itervalues
from six.moves import xrange as range

from cellranger.library_constants import GENE_EXPRESSION_LIBRARY_TYPE, ATACSEQ_LIBRARY_TYPE
from atac.cell_calling_helpers.exclusions import LOW_TARGETING, GEL_BEAD_DOUBLET

MAX_CELLS = 80000    # Overwrite for overloaded 10x multiome samples

# these are ATAC barcode exclusions that we consider prior to joint cell calling
JOINT_EXCLUSION_REASONS = [GEL_BEAD_DOUBLET, LOW_TARGETING]


class CellCallerException(Exception):
    """ Exception class for cell caller errors """


def parse_atac_exclusions(excluded_barcodes, genomes):
    """
    Outputs a dictionary of barcodes excluded due to gel bead doublets or excess whitelist
    contamination keyed by genome.
    """
    exclusions = {}

    for genome in genomes:
        # for a single genome, the excluded barcodes use a key "", while for multiple genomes
        # the key is the genome name
        if len(genomes) == 1:
            gdata = next(itervalues(excluded_barcodes))
        else:
            gdata = excluded_barcodes[genome]
        bcs = set()
        for ex_bc, (reason, _) in iteritems(gdata):
            if reason in JOINT_EXCLUSION_REASONS:
                bcs.add(ensure_binary(ex_bc))
        exclusions[genome] = bcs

    return exclusions


def construct_counts_dataframe(rna_matrix, atac_matrix, genomes, excluded_barcodes=None):
    """
    For each genome construct a pandas dataframe that is the input to the cell caller.
    excluded_barcodes is a dictionary {genome: set(excluded barcodes)}. The result is a dict
    {genome: dataframe} where the dataframe.index = barcode and the atac and rna counts are
    columns names 'atac_count' and 'rna_count' respectively. The barcode exclusions are in the
    'excluded' column.
    """
    df_per_genome = {}
    for genome in genomes:
        rna_count = rna_matrix.get_counts_per_barcode_for_genome(
            genome, GENE_EXPRESSION_LIBRARY_TYPE
        )
        non_zero_counts = rna_count > 0
        rna_df = pd.DataFrame()
        rna_df["barcode"] = rna_matrix.bcs[non_zero_counts]
        rna_df["rna_count"] = rna_count[non_zero_counts]

        # Note: we allow barcodes with zero atac counts here since some of these could appear as
        # excluded barcodes below
        atac_df = pd.DataFrame()
        atac_df["barcode"] = atac_matrix.bcs
        atac_df["atac_count"] = atac_matrix.get_counts_per_barcode_for_genome(
            genome, ATACSEQ_LIBRARY_TYPE
        )

        df = pd.merge(rna_df, atac_df, on="barcode", how="outer")
        df.fillna(value=0, inplace=True)
        df.set_index("barcode", inplace=True)

        # after merge the counts become floats
        df["atac_count"] = df["atac_count"].astype(np.int32)
        df["rna_count"] = df["rna_count"].astype(np.int32)

        df["excluded"] = False
        if excluded_barcodes:
            ebcs = excluded_barcodes.get(genome)
            if ebcs:
                df.loc[ebcs, "excluded"] = True
        df_per_genome[genome] = df
    return df_per_genome


def filter_by_min_counts(data, min_atac, min_gex):
    """Outputs a mask of barcodees whose ATAC AND GEX counts
    are above the input minimum for each assay
    """
    return np.logical_and(data["atac_count"] >= min_atac, data["rna_count"] >= min_gex).astype(
        "int"
    )


def simpler_ordmag(x):
    """
    Do an initial estimate of ordmag threshold using 99% of all the barcodes. Of the barcodes
    greater than the threshold discard the top 1%. Repeat ordmag on remaining barcodes to define
    the final threshold.
    """
    assert len(x) > 0, "Must specify a non-zero length vector for running ordmag"
    p99_1 = np.percentile(x, 99)
    ordmag = x > p99_1 / 10.0
    # this should never happen
    if ordmag.sum() == 0:
        return 0
    top_1_ordmag = np.percentile(x[ordmag], 99)
    ordmag2 = x < top_1_ordmag
    # use the previous threshold
    if ordmag2.sum() == 0:
        return p99_1 / 10.0
    p99_2 = np.percentile(x[ordmag2], 99)
    return p99_2 / 10.0


def log_transform(x, alpha=0.1):
    return np.log(x + alpha)


def transfer_labels_to_duplicates(is_duplicate, values):
    """
    Given a set of values and a boolean array marking duplicates, copy the values from a
    non-duplicate entry to all duplicate entries that immediately follow.
    """
    assert is_duplicate.shape == values.shape
    non_dup = ~is_duplicate
    values_dedup = values[non_dup]
    dedup_ptr = np.cumsum(non_dup) - 1
    return values_dedup[dedup_ptr]


def cap_cells(data, max_cells):
    """
    If more than max_cells have been called, restrict to the top max_cells barcodes based on
    atac_count + rna_count.
    """
    if data["is_cell"].sum() > max_cells:
        joint_count = data["atac_count"] + data["rna_count"]
        joint_count = joint_count.sort_values(ascending=False).index
        sorted_data = data.loc[joint_count, ["is_cell"]]
        discard_cells = sorted_data[sorted_data["is_cell"] == 1].index[max_cells:]
        data.loc[discard_cells, "is_cell"] = 0


def define_keep_filter(data, count_threshold):
    """ Which barcodes are considered for cell calling """
    if data["excluded"].sum() == data.shape[0]:
        raise CellCallerException(
            "All barcodes have been excluded from cell calling due to low-targeting and/or being "
            "gel bead doublets. This exclusion is based on the ATAC data alone and could be "
            "caused by a failure to identify peaks in the ATAC data or an error in sample/library "
            "preparation."
        )

    rna_filter = data["rna_count"] > count_threshold
    atac_filter = data["atac_count"] > count_threshold
    low_count_msg = (
        "All barcodes have <= {} counts for {}. Cell calling cannot be performed with low counts. "
        "This could be caused by a failure in library preparation or due to very low sequencing "
        "depth for the {}."
    )
    if atac_filter.sum() == 0:
        raise CellCallerException(low_count_msg.format(count_threshold, "ATAC", "ATAC library"))

    if rna_filter.sum() == 0:
        raise CellCallerException(low_count_msg.format(count_threshold, "GEX", "GEX library"))

    if (atac_filter & rna_filter).sum() == 0:
        raise CellCallerException(
            "Barcodes must have > {} counts for both ATAC and GEX modalities. This could be "
            "caused by a failure in library preparation or very low sequencing depth or "
            "a mis-pairing of ATAC library FASTQs from one sample with GEX library FASTQs of a "
            "different one".format(count_threshold)
        )
    keep = (~data["excluded"]) & (~data["dup"]) & rna_filter & atac_filter

    if keep.sum() == 0:
        raise CellCallerException(
            "No barcodes left after excluding barcodes based on the ATAC data as low targeting or "
            "gel bead doublets and requiring at least {} count per barcode for both ATAC and GEX "
            "modalities. This could be caused by very low sequencing depth, or an error during "
            "sample/library preparation, or a failure to identify peaks in the ATAC data, or "
            "mis-pairing of ATAC library FASTQs from one sample with GEX library FASTQs of a "
            "different one.".format(count_threshold)
        )
    return keep


def call_cells(data, force_cells, count_threshold=1, max_cells=MAX_CELLS):
    """
    Inputs:
        data: dataframe with columns 'atac_counts' (int), 'rna_counts' (int), 'excluded' (bool).
        force_cells: {"atac": min_atac_count, "gex": min_gex_count} for forcing cells. Note that
                     barcodes will be excluded before forcing. Can be null.
        count_threshold: any barcode with < this threshold counts in ATAC or GEX is guaranteed to
                         not be a cell.
        max_cells: do not call more than max_cells cells.
    Outputs:
        parameters: dict of cell calling metrics

    Note: that the dataframe supplied is modified in place to add new columns `excluded`, `is_cell`,
          and `dup`.

    1. Only consider de-duplicated barcodes, barcodes that are not excluded by atac and with counts
       above the count_threshold in both assays
    2. Determine initial centroids of cell and non-cell clouds using ord-mag
    3. Run K-means with K=2 to refine

    if force_cells is specified, override this behavior using the min atac and rna counts specified.
    """

    data["is_cell"] = 0

    # sort values to put duplicates next to each other
    data.sort_values(by=["atac_count", "rna_count", "excluded"], inplace=True)
    data["dup"] = data.duplicated(keep="first")

    parameters = {}
    for assay in ["rna", "atac"]:
        parameters[assay] = {
            "ordmag_threshold": np.nan,
            "centroid_initial": np.nan,
            "centroid_final": np.nan,
            "force_cells_min_count": np.nan,
        }

    # force cells if specified, but do not select excluded barcodes
    if force_cells:
        data["is_cell"] = filter_by_min_counts(data, force_cells["atac"], force_cells["gex"])
        parameters["rna"]["force_cells_min_count"] = force_cells["gex"]
        parameters["atac"]["force_cells_min_count"] = force_cells["atac"]
        data.loc[data["excluded"], "is_cell"] = 0
    else:
        # only consider non-excluded and non-duplicated barcodes with > count_threshold counts
        keep = define_keep_filter(data, count_threshold)

        # 1: do ordmag calling and identify centroids
        threshold = {}
        for assay in ["atac", "rna"]:
            threshold[assay] = simpler_ordmag(data.loc[keep, "{}_count".format(assay)])
            parameters[assay]["ordmag_threshold"] = threshold[assay]
        ordmag_filter = (
            (data["atac_count"] >= threshold["atac"]) & (data["rna_count"] >= threshold["rna"])
        ).astype(np.int8)

        clf = NearestCentroid()
        points = log_transform(data.loc[keep, ["atac_count", "rna_count"]].values)
        clsfn = np.array(ordmag_filter[keep])

        # If ordmag finds that all barcodes are to be kept, this is an edge case, and call all
        # barcodes as cells. The number of cells will be capped by `cap_cells()` below. If we
        # call everything a cell, the downstream analyses will at least tell us something about
        # the high count barcodes in the data.
        if len(set(clsfn)) == 1:
            data.loc[keep, "is_cell"] = 1
        else:
            _ = clf.fit(points, clsfn)
            for i, assay in enumerate(["atac", "rna"]):
                parameters[assay]["centroid_initial"] = np.exp(clf.centroids_[:, i]).tolist()

            # 2: k-means initialized with centroids from above
            kmeans_labels = KMeans(n_clusters=2, init=clf.centroids_, n_init=1).fit_predict(points)

            # if kmeans puts all the points in one cluster then call everything a cell
            if len(set(kmeans_labels)) == 1:
                data.loc[keep, "is_cell"] = 1
            else:
                # make sure labels didn't get switched
                mean0 = points[kmeans_labels == 0, :].mean()
                mean1 = points[kmeans_labels == 1, :].mean()
                if mean1 < mean0:
                    kmeans_labels = 1 - kmeans_labels
                data.loc[keep, "is_cell"] = kmeans_labels

                clf.fit(points, kmeans_labels)
                for i, assay in enumerate(["atac", "rna"]):
                    parameters[assay]["centroid_final"] = np.exp(clf.centroids_[:, i]).tolist()

        # 3: transfer labels to duplicates
        data["is_cell"] = transfer_labels_to_duplicates(data["dup"].values, data["is_cell"].values)

    # Ensure that every called cell has at least count_threshold UMIs and fragments
    zero_signal = np.logical_or(
        data["rna_count"] < count_threshold, data["atac_count"] < count_threshold
    )
    data.loc[zero_signal, "is_cell"] = 0

    # Make sure we don't call more than max_cells cells
    cap_cells(data, max_cells)

    data["dup"] = data["dup"].astype(np.int8)

    return parameters


def simulate_data(spec):
    """
    Simulate cell calling data based on a spec. An example spec is as follows:
    [
        {
            "name": "contam",
            "points": 10000,
            "atac": {"mean": 1, "spread": 1},
            "rna": {"mean": 1, "spread": 1},
        },
        {
            "name": "cells",
            "points": 100,
            "atac": {"mean": 10000, "spread": 1},
            "rna": {"mean": 10000, "spread": 1},
        },
        {
            "name": "pop",
            "points": 10,
            "atac": {"mean": 400, "spread": 1},
            "rna": {"mean": 400, "spread": 1},
        },
    ]
    Each element specifies a population with a mean and spread for both assays and how many
    points to generate.
    """

    def make_counts(mean, spread, num_points):
        return np.random.lognormal(np.log(mean), spread, size=num_points).astype(np.int32)

    num_barcodes = sum(x["points"] for x in spec)
    atac = np.zeros(num_barcodes, dtype=np.int32)
    rna = np.zeros(num_barcodes, dtype=np.int32)

    truth = []
    index = 0
    for btype in spec:
        points = btype["points"]
        atac[index : index + points] = make_counts(
            btype["atac"]["mean"], btype["atac"]["spread"], points
        )
        rna[index : index + points] = make_counts(
            btype["rna"]["mean"], btype["rna"]["spread"], points
        )
        index += points
        truth.extend([btype["name"] for _ in range(points)])

    data = pd.DataFrame()
    data["atac_count"] = atac
    data["rna_count"] = rna
    data["truth"] = truth
    data["excluded"] = False
    return data
