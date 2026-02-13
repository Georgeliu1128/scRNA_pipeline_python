import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
import os
import logging
from datetime import datetime


# ---------- LOGGING SETUP ----------
os.makedirs("logs", exist_ok=True)

logfile = f"logs/pipeline_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    filename=logfile,
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

logger = logging.getLogger(__name__)
logger.info("Pipeline started")




sc.settings.autosave = True
sc.settings.autoshow = False
sc.settings.figdir = "figures"
sc.settings.file_format_figs = "png"
sc.settings.verbosity = 2
sc.set_figure_params(figsize=(5, 5), dpi=100)



def savefig(name, out_prefix="run"):
    os.makedirs("figures", exist_ok=True)
    plt.savefig(f"figures/{out_prefix}_{name}.png", bbox_inches="tight")
    plt.close()



# ---------- IO / CONCAT ----------

def load_and_concat(h5_paths, sample_names):
    adatas = []
    for path, name in zip(h5_paths, sample_names):
        adata = sc.read_10x_h5(path)
        adata.var_names_make_unique()     # <-- FIX
        adata.obs_names_make_unique()     # <-- FIX (rare but safe)
        adata.obs["sample"] = name
        adatas.append(adata)
    return ad.concat(adatas, label="batch", keys=sample_names, index_unique="-")



# ---------- QC MODULE ----------

def run_qc(
    adata,
    mt_prefix="MT-",
    mt_thresh=15,
    min_genes=300,
    max_genes=6000,
    min_cells=3,
    out_prefix="qc"
):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(mt_prefix)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # Scanpy autosaves these automatically
    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        groupby="batch",
        multi_panel=True
    )

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", color="batch")
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="batch")

    # Filtering
    adata = adata[adata.obs["pct_counts_mt"] < mt_thresh, :]
    adata = adata[adata.obs["n_genes_by_counts"] > min_genes, :]
    adata = adata[adata.obs["n_genes_by_counts"] < max_genes, :]

    sc.pp.filter_genes(adata, min_cells=min_cells)
    return adata



# ---------- NORMALIZATION + HVG MODULE ----------

def normalize_and_hvg(
    adata,
    n_top_genes=3000,
    batch_key="batch",
    out_prefix="hvg"
):
    """Library size normalize, log1p, HVG selection."""
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor="seurat_v3",
        batch_key=batch_key
    )

    sc.pl.highly_variable_genes(adata)

    adata = adata[:, adata.var["highly_variable"]].copy()
    return adata


# ---------- BATCH CORRECTION BACKENDS ----------

def run_pca(adata, n_pcs=50):
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver="arpack", n_comps=n_pcs)
    return adata


def batch_correct_none(adata, n_pcs=50):
    """No batch correction, just PCA."""
    adata = run_pca(adata, n_pcs=n_pcs)
    return adata, "X_pca"


def batch_correct_harmony(adata, batch_key="batch", n_pcs=50):
    import harmonypy as hm
    import numpy as np

    # Run PCA
    adata = run_pca(adata, n_pcs=n_pcs)

    # Extract PCA matrix (cells × pcs)
    Z = adata.obsm["X_pca"]

    # --- CRITICAL METADATA CLEANING ---
    meta = adata.obs[[batch_key]].copy()

    # Force batch column to string categories
    meta[batch_key] = meta[batch_key].astype(str)

    # Reset index completely (Harmony hates named or non-string indices)
    meta.index = meta.index.astype(str)
    meta.index.name = None

    # Run Harmony
    ho = hm.run_harmony(Z, meta, batch_key)

    # Extract corrected matrix
    Z_corr = np.asarray(ho.Z_corr)

    # If Harmony returned (pcs × cells), transpose it
    if Z_corr.shape == (n_pcs, adata.n_obs):
        Z_corr = Z_corr.T

    # Final shape validation
    if Z_corr.shape != Z.shape:
        raise ValueError(
            f"Harmony output shape mismatch: expected {Z.shape}, got {Z_corr.shape}"
        )

    adata.obsm["X_harmony"] = Z_corr
    return adata, "X_harmony"



def batch_correct_bbknn(adata, batch_key="batch", n_pcs=50):
    """BBKNN batch correction (neighbors in batch-aware space)."""
    import bbknn

    adata = run_pca(adata, n_pcs=n_pcs)
    bbknn.bbknn(adata, batch_key=batch_key, n_pcs=n_pcs)
    # BBKNN writes neighbors directly; use X_pca for embedding
    return adata, "X_pca"


# ---------- NEIGHBORS / UMAP / CLUSTERING MODULE ----------

def neighbors_umap_cluster(
    adata,
    rep_key,
    n_neighbors=15,
    n_pcs=40,
    resolution=0.5,
    out_prefix="embed"
):
    """Build graph, UMAP, Leiden clustering from a given representation."""
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, use_rep=rep_key)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=resolution)

    # UMAP plots
    sc.pl.umap(adata, color=["batch"])
    sc.pl.umap(adata, color=["sample"])
    sc.pl.umap(adata, color=["leiden"], legend_loc="on data")

    return adata


# ---------- DE / MARKERS MODULE ----------

def run_markers(
    adata,
    groupby="leiden",
    method="wilcoxon",
    out_prefix="markers"
):
    # --- RUN DIFFERENTIAL EXPRESSION FIRST ---
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method
    )

    # --- THEN PLOT ---
    sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
    sc.pl.rank_genes_groups_heatmap(adata, n_genes=5, groupby=groupby, dendrogram=False)

    return adata


# ---------- TOP-LEVEL PIPELINE WRAPPER ----------

def run_full_pipeline(
    h5_paths,
    sample_names,
    batch_backend="none",  # "none", "harmony", "bbknn"
    out_prefix="run"
):
    adata = load_and_concat(h5_paths, sample_names)
    adata = run_qc(adata, out_prefix=f"{out_prefix}_qc")
    adata = normalize_and_hvg(adata, out_prefix=f"{out_prefix}_hvg")

    if batch_backend == "none":
        adata, rep = batch_correct_none(adata)
    elif batch_backend == "harmony":
        adata, rep = batch_correct_harmony(adata)
    elif batch_backend == "bbknn":
        adata, rep = batch_correct_bbknn(adata)
    else:
        raise ValueError(f"Unknown batch backend: {batch_backend}")

    adata = neighbors_umap_cluster(
        adata,
        rep_key=rep,
        out_prefix=f"{out_prefix}_embed"
    )
    adata = run_markers(adata, out_prefix=f"{out_prefix}_markers")

    adata.write_h5ad(f"{out_prefix}_processed.h5ad")
    return adata
