# Python wrapper functions for the Spatial transcriptomics pipeline

from typing import Literal, Union, Optional
import scanpy as sc
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import squidpy as sq
# import decoupler as dc
from pathlib import Path
from anndata import AnnData
from enum import Enum, unique

# TODO: The functions need to be documented

@unique
class Protocol(Enum):
    FF = "FF"
    FFPE = "FFPE"
    
@unique
class Organism(Enum):
    mouse = "mouse"
    human = "human"
    rat = "rat"
    

class Analyze():
    def __init__(self, protocol: Protocol, organism: Organism):
        self.protocol = protocol
        self.organism = organism

        if self.organism.value == Organism.human.value:
            self.mito_prefix = "MT-"
        elif self.organism.value == Organism.mouse.value:
            self.mito_prefix = "mt-"
        elif self.organism.value == Organism.rat.value:
            ## Modified Alberto
            self.mito_prefix = "Mt-"
        else:
            self.mito_prefix = ""
        
    def __str__(self):
        return(f"Protocol: {self.protocol}, Organism: {self.organism}, Mitochondrial Gene Prefix: {self.mito_prefix}")


##########################
# General functions to import data
# #########################


def generate_adata_objects(
    path: Union[str, Path],
    samples_names: list[str],
    metadata: pd.DataFrame(),
    analyze_params: Analyze,
    count_file: str = "filtered_feature_bc_matrix.h5",
    load_images: Optional[bool] = True,
) -> list[AnnData]:

    adata_objects = []

    for current_sample_name in samples_names:

        current_adata = sq.read.visium(
            path=os.path.join(path, current_sample_name),
            counts_file=count_file,
            load_images=load_images
            )
        current_adata.var_names_make_unique()
        current_adata.obs["Sample_ID"] = current_sample_name

        for column in metadata.columns:
            current_adata.obs[column] = metadata.loc[current_sample_name, column]

        # It may be that in the new probe sets versions for the FFPE protocol there are some 
        # mito genes, as we have seen in singlecell Flex experiements. I would therefore
        # remove the "if" for now for simplicity. To monitor potential differences between FF and
        # FFPE that may be worthy to include. 
        # if analyze_params.protocol.value == Protocol.FF.value:
        #    print("Fresh frozen samples: Detecting mitochondrial genes")
        current_adata.var["mt"] = current_adata.var_names.str.startswith(analyze_params.mito_prefix)
        sc.pp.calculate_qc_metrics(current_adata, qc_vars=["mt"], inplace=True)

        adata_objects.append(current_adata)

    return adata_objects


##########################
# QC related functions
# #########################


def get_global_QCmetrics(
    path: str, samples_names: list[str], metrics_names: str = "metrics_summary.csv"
) -> pd.DataFrame():

    global_metrics = pd.DataFrame()
    for current_sample in samples_names:
        current_metrics = pd.read_csv(
            os.path.join(path, current_sample, metrics_names)
        )
        global_metrics = pd.concat([global_metrics, current_metrics])

    return global_metrics.set_index("Sample ID")


# We can simplify depending if the batch number comes in the name in a standard format
def get_barplot_qc(
    qc_df: pd.DataFrame(),
    color_by: list[str],
    variables_to_plot: np.ndarray,
    plots_row: Optional[int] = 3,
):

    for current_variable in variables_to_plot:

        current_df = qc_df[[current_variable]]
        current_df = current_df.assign(color_by=color_by)

        sns.set_theme(style="darkgrid")
        sns.barplot(
            data=current_df,
            x=current_variable,
            y=current_df.index,
            hue=current_df["color_by"],
            dodge=False,
        )
        # g.set_xticklabels(g.get_xticklabels(), rotation=30)
        plt.show()
        plt.close()


def perform_qc_analysis(
    list_adata_filter: list[AnnData],
    list_adata_raw: list[AnnData],
    color_map="OrRd",
    sample_id="Sample_ID",
    condition_name="Condition",
    batch_name="Batch"
    
):

    for a in range(len(list_adata_filter)):

        id = "".join(set(list_adata_raw[a].obs[sample_id]))
        condition = "".join(set(list_adata_raw[a].obs[condition_name]))
        batch = "".join(set(list_adata_raw[a].obs[batch_name]))
        title = f"{sample_id}: {id}, {condition_name}: {condition}, {batch_name}: {batch}"
        sc.pl.spatial(adata=list_adata_raw[a], img_key="hires", title=title, show=False)

        print(f"All spots: {title}")
        sc.pl.spatial(
            list_adata_raw[a],
            img_key="hires",
            color=["total_counts", "n_genes_by_counts"],
            color_map=color_map,
        )

        print(f"Tissue covered spots: {title}")
        sc.pl.spatial(
            list_adata_filter[a],
            img_key="hires",
            color=["total_counts", "n_genes_by_counts"],
            color_map=color_map,
        )

        fig, axs = plt.subplots(1, 2, figsize=(15, 5))
        sns.violinplot(
            data=list_adata_raw[a].obs,
            x="in_tissue",
            y="total_counts",
            inner="quart",
            linewidth=1,
            ax=axs[0],
        )
        sns.stripplot(
            y="total_counts",
            x="in_tissue",
            data=list_adata_raw[a].obs,
            color="black",
            edgecolor="gray",
            size=2.5,
            ax=axs[0],
        )

        sns.violinplot(
            data=list_adata_raw[a].obs,
            x="in_tissue",
            y="n_genes_by_counts",
            inner="quart",
            linewidth=1,
            ax=axs[1],
        )
        sns.stripplot(
            y="n_genes_by_counts",
            x="in_tissue",
            data=list_adata_raw[a].obs,
            color="black",
            edgecolor="gray",
            size=2.5,
            ax=axs[1],
        )
        plt.show(fig)
        plt.close(fig)

        fig, axs = plt.subplots(1, 3, figsize=(22, 5))
        sns.scatterplot(
            data=list_adata_filter[a].obs,
            x="total_counts",
            y="n_genes_by_counts",
            ax=axs[0],
        )
        sns.violinplot(
            data=list_adata_filter[a].obs,
            x="Sample_ID",
            y="total_counts",
            inner="quart",
            linewidth=1,
            ax=axs[1],
        )
        sns.stripplot(
            y="total_counts",
            x="Sample_ID",
            data=list_adata_filter[a].obs,
            color="black",
            edgecolor="gray",
            size=2.5,
            ax=axs[1],
        )

        sns.violinplot(
            data=list_adata_filter[a].obs,
            x="Sample_ID",
            y="n_genes_by_counts",
            inner="quart",
            linewidth=1,
            ax=axs[2],
        )
        sns.stripplot(
            y="n_genes_by_counts",
            x="Sample_ID",
            data=list_adata_filter[a].obs,
            color="black",
            edgecolor="gray",
            size=2.5,
            ax=axs[2],
        )
        plt.show(fig)
        plt.close(fig)

        # TODO: Check for the FFPE Implementation
        print(f"Mithocondrial genes: {title}")
        sc.pl.spatial(
            list_adata_filter[a],
            img_key="hires",
            color=["total_counts_mt", "pct_counts_mt"],
            color_map=color_map,
        )

        fig, axs = plt.subplots(1, 3, figsize=(22, 5))
        sns.scatterplot(
            data=list_adata_filter[a].obs,
            x="total_counts",
            y="total_counts_mt",
            ax=axs[0],
        )

        sns.violinplot(
            data=list_adata_filter[a].obs,
            x="Sample_ID",
            y="total_counts_mt",
            inner="quart",
            linewidth=1,
            ax=axs[1],
        )
        sns.stripplot(
            y="total_counts_mt",
            x="Sample_ID",
            data=list_adata_filter[a].obs,
            color="black",
            edgecolor="gray",
            size=2.5,
            ax=axs[1],
        )

        sns.violinplot(
            data=list_adata_filter[a].obs,
            x="Sample_ID",
            y="pct_counts_mt",
            inner="quart",
            linewidth=1,
            ax=axs[2],
        )
        sns.stripplot(
            y="pct_counts_mt",
            x="Sample_ID",
            data=list_adata_filter[a].obs,
            color="black",
            edgecolor="gray",
            size=2.5,
            ax=axs[2],
        )
        plt.show(fig)
        plt.close(fig)


# By default, I am entering very "relaxed" thresholds
def qc_filtering(
    list_adatas: list[AnnData],
    min_counts: Optional[int] = 1000,
    max_counts: Optional[int] = 40000,
    mt_pct_content: Optional[int] = 20,
    min_cells: Optional[int] = 5,
) -> list[AnnData]:

    adata_filtered_objects: list[AnnData] = []

    for current_adata in list_adatas:

        current_sample = np.asarray(current_adata.obs["Sample_ID"].unique())
        print(current_sample)

        print(f"# Spots before filter: {current_adata.n_obs}")
        print(f"# Genes before filter: {current_adata.n_vars}")

        sc.pp.filter_cells(current_adata, min_counts=min_counts, inplace=True)
        sc.pp.filter_cells(current_adata, max_counts=max_counts, inplace=True)
        current_adata = current_adata[
            current_adata.obs["pct_counts_mt"] < mt_pct_content
        ]
        sc.pp.filter_genes(current_adata, min_cells=min_cells, inplace=True)

        print(f"# Spots after filter: {current_adata.n_obs}")
        print(f"# Genes before filter: {current_adata.n_vars}")

        adata_filtered_objects.append(current_adata)

    return adata_filtered_objects


##########################
# Normalization, manifold embedding, clustering and Marker genes
# #########################

# TODO: more parameters to be added.
def norm_hvg(
    list_adatas: list[AnnData],
    flavor: Literal["seurat", "cell_ranger", "seurat_v3"] = "seurat",
    n_top_genes: Optional[int] = None,
) -> list[AnnData]:

    adata_objects: list[AnnData] = []

    for adata in list_adatas:
        sc.pp.normalize_total(adata, inplace=True)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=n_top_genes)
        adata_objects.append(adata)

    return adata_objects


def cluster_umap(
    list_adatas: list[AnnData],
    n_comps: Optional[int] = None,
    use_highly_variable: Optional[bool] = None,
    n_neighbors: Optional[int] = 15,
    n_pcs: Optional[int] = None,
    min_dist: Optional[float] = 0.5,
    spread: Optional[float] = 1.0,
    n_components: Optional[int] = 2,
    key_added: str = "clusters",
    resolution: float = 0.75,
) -> list[AnnData]:

    adata_objects: list[AnnData] = []

    for current_adata in list_adatas:

        current_sample = np.asarray(current_adata.obs["Sample_ID"].unique())
        print(current_sample)

        sc.pp.pca(
            current_adata, n_comps=n_comps, use_highly_variable=use_highly_variable
        )
        sc.pp.neighbors(current_adata, n_neighbors=n_neighbors, n_pcs=n_pcs)
        sc.tl.umap(
            current_adata, min_dist=min_dist, spread=spread, n_components=n_components
        )
        sc.tl.leiden(current_adata, resolution=resolution, key_added=key_added)

        plt.rcParams["figure.figsize"] = (4, 4)
        sc.pl.umap(
            current_adata,
            color=["total_counts", "n_genes_by_counts", "clusters"],
            wspace=0.4,
        )

        plt.rcParams["figure.figsize"] = (6, 6)
        sc.pl.spatial(current_adata, img_key="hires", color="clusters", size=1.5)

        adata_objects.append(current_adata)

    return adata_objects


def get_markers_clusters(
    list_adatas,
    group_by="clusters",
    groups_1="all",
    reference="rest",
    method="wilcoxon",
    corr_method="benjamini-hochberg",
    groups_2=None,
    n_genes=5,
    genes_visualize=20,
) -> list[AnnData]:

    adata_objects: list[AnnData] = []

    for current_adata in list_adatas:

        current_sample = np.asarray(current_adata.obs["Sample_ID"].unique())
        print(current_sample)

        sc.tl.rank_genes_groups(
            current_adata,
            group_by,
            groups=groups_1,
            reference=reference,
            method=method,
            corr_method=corr_method,
        )
        plt.rcParams["figure.figsize"] = (6, 6)
        sc.pl.rank_genes_groups_heatmap(
            current_adata, groupby=group_by, groups=groups_2, n_genes=n_genes
        )

        # plt.rcParams["figure.figsize"] = (6, 6)
        # sc.pl.spatial(current_adata, img_key="hires", color=["clusters"], size = 1.5)

        group = list(current_adata.uns["rank_genes_groups"]["names"].dtype.names)
        colnames = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]
        d = [
            pd.DataFrame(current_adata.uns["rank_genes_groups"][c])[group]
            for c in colnames
        ]
        d = pd.concat(d, axis=1, names=[None, "group"], keys=colnames)
        d = d.stack(level=1).reset_index()
        d["group"] = pd.Categorical(d["group"], categories=group)
        d = d.drop(columns="level_0")
        d = d.sort_values(["pvals_adj"], ascending=True)

        print(d.head(genes_visualize))

        adata_objects.append(current_adata)

    return adata_objects


##########################
# Spatially variable genes
# #########################


# max_neighs Only used for Septal
def get_sp_variable_genes(
    list_adatas,
    method="moran",
    min_number_spots=100,
    number_hvg=100,
    n_perms=100,
    n_jobs=1,
    genes_visualize=10,
    max_neighs=6,
) -> list[AnnData]:

    adata_objects: list[AnnData] = []

    for current_adata in list_adatas:

        current_sample = np.asarray(current_adata.obs["Sample_ID"].unique())
        print(current_sample)

        genes = current_adata.var_names[
            (current_adata.var.n_cells > min_number_spots)
            & current_adata.var.highly_variable
        ][0:number_hvg]
        sq.gr.spatial_neighbors(current_adata)
        genes = current_adata.var_names[
            (current_adata.var.n_cells > min_number_spots)
            & current_adata.var.highly_variable
        ][0:number_hvg]

        if method == "sepal":
            sq.gr.sepal(
                current_adata, max_neighs=max_neighs, genes=genes, n_jobs=n_jobs
            )
            print(current_adata.uns["sepal_score"].head(genes_visualize))
        elif method == "moran":
            sq.gr.spatial_autocorr(
                current_adata, mode="moran", genes=genes, n_perms=n_perms, n_jobs=n_jobs
            )
            print(current_sample)
            print(current_adata.uns["moranI"].head(genes_visualize))
        else:
            print("Unknown method")

        adata_objects.append(current_adata)

    return adata_objects


##########################
# Footprint-based methods
# #########################


def get_pathway_activity(
    list_adatas,
    organism="human",
    top_genes=500,
    verbose=False,
    groupby="clusters",
    vmin=-2,
    vmax=2,
    cmap="coolwarm",
    use_raw=False,
) -> list[AnnData]:

    adata_objects: list[AnnData] = []
    model = dc.get_progeny(organism=organism, top=top_genes)

    for a in range(len(list_adatas)):

        current_sample = np.asarray(list_adatas[a].obs["Sample_ID"].unique())
        print(current_sample)

        current_adata = list_adatas[a]

        dc.run_mlm(
            mat=current_adata,
            net=model,
            source="source",
            target="target",
            weight="weight",
            verbose=verbose,
            use_raw=use_raw,
        )

        current_adata.obsm["progeny_mlm_estimate"] = current_adata.obsm[
            "mlm_estimate"
        ].copy()
        current_adata.obsm["progeny_mlm_pvals"] = current_adata.obsm["mlm_pvals"].copy()

        acts = dc.get_acts(current_adata, obsm_key="progeny_mlm_estimate")

        mean_acts = dc.summarize_acts(acts, groupby=groupby, min_std=0)
        print(mean_acts)

        sns.clustermap(
            mean_acts, xticklabels=mean_acts.columns, vmin=vmin, vmax=vmax, cmap=cmap
        )
        plt.show()

        adata_objects.append(current_adata)

    return adata_objects


def get_TF_activity(
    list_adatas,
    organism="human",
    levels=["A", "B", "C"],
    min_n=5,
    min_std=0.75,
    verbose=False,
    groupby="clusters",
    vmin=-2,
    vmax=2,
    cmap="coolwarm",
    use_raw=False,
) -> list[AnnData]:

    adata_objects: list[AnnData] = []
    net = dc.get_dorothea(organism=organism, levels=["A", "B", "C"])

    for a in range(len(list_adatas)):

        current_sample = np.asarray(list_adatas[a].obs["Sample_ID"].unique())
        print(current_sample)

        current_adata = list_adatas[a]

        dc.run_mlm(
            mat=current_adata,
            net=net,
            min_n=min_n,
            source="source",
            target="target",
            weight="weight",
            verbose=verbose,
            use_raw=use_raw,
        )

        current_adata.obsm["dorothea_mlm_estimate"] = current_adata.obsm[
            "mlm_estimate"
        ].copy()
        current_adata.obsm["dorothea_mlm_pvals"] = current_adata.obsm[
            "mlm_pvals"
        ].copy()

        acts = dc.get_acts(current_adata, obsm_key="dorothea_mlm_estimate")

        mean_acts = dc.summarize_acts(acts, groupby=groupby, min_std=min_std)
        print(mean_acts)

        sns.clustermap(
            mean_acts, xticklabels=mean_acts.columns, vmin=vmin, vmax=vmax, cmap=cmap
        )
        plt.show()

        adata_objects.append(current_adata)

    return adata_objects
