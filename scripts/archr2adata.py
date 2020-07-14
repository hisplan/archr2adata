#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import scanpy as sc
import scipy.io

def make_adata(path_data):

    # load
    barcodes = pd.read_csv(
        os.path.join(path_data, "peak_counts/cells.csv.gz"),
        skiprows=[0], header=None
    )[1]

    peaks = pd.read_csv(
        os.path.join(path_data, "peak_counts/peaks.csv.gz"),
        index_col=0
    )

    counts = scipy.io.mmread(
        os.path.join(path_data, "peak_counts/counts.mtx.gz")
    )

    # construct AnnData
    adata = sc.AnnData(counts.T)

    adata.obs_names = barcodes
    adata.var_names = peaks.index

    for col in ['distToGeneStart', 'nearestGene', 'peakType', 'distToTSS', 'nearestTSS','replicateScoreQuantile']:
        adata.var[col] = peaks[col]

    # SVD
    adata.obsm["X_svd"] = pd.read_csv(
        os.path.join(path_data, "svd.csv"),
        index_col=0
    ).values

    # UMAP
    adata.obsm["umap"] = pd.read_csv(
        os.path.join(path_data, "umap.csv"),
        index_col=0
    ).values

    # metadata
    df_meta = pd.read_csv(
        os.path.join(path_data, "cell_metadata.csv"),
        index_col=0
    )

    for col in df_meta.columns:
        adata.obs[col] = df_meta[col].values

    # gene scores
    df_gene_scores = pd.read_csv(
        os.path.join(path_data, "gene_scores.csv.gz"),
        index_col=0
    ).T

    adata.obsm['GeneScores'] = df_gene_scores.values
    adata.uns['GeneScoresColumns'] = df_gene_scores.columns.values

    # a writer for type <class 'scipy.sparse.coo.coo_matrix'> has not been implemented yet.
    from scipy.sparse import csr_matrix

    adata.X = csr_matrix(adata.X)

    adata.write(
        os.path.join(path_data, "preprocessed.h5ad")
    )


def parse_arguments():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--data",
        action="store",
        dest="path_data",
        help="path to data files (*.csv, *.mtx, ...)",
        required=True
    )

    # parse arguments
    params = parser.parse_args()

    return params


if __name__ == "__main__":

    params = parse_arguments()

    make_adata(params.path_data)
