from bson import decode_file_iter
import pandas as pd


def bson_to_df(path: str) -> pd.DataFrame:
    docs = []
    with open(path, "rb") as f:
        for d in decode_file_iter(f):
            docs.append(d)
    return pd.json_normalize(docs, sep=".")


def has_product(x) -> bool:
    return isinstance(x, list) and len(x) > 0


def is_list(x) -> bool:
    return isinstance(x, list)


# ----------------------------
# Load data
# ----------------------------
gene_df = bson_to_df("../data/regulonDB_14/regulondbdatamarts/geneDatamart.bson")
regulon_df = bson_to_df("../data/regulonDB_14/regulondbdatamarts/regulonDatamart.bson")
tu_df = bson_to_df("../data/regulonDB_14/regulondbht/transcriptionUnit.bson")


# ----------------------------
# Gene universe (protein-coding proxy)
# ----------------------------
gene_universe = set(gene_df["gene._id"].astype(str))
protein_genes = set(gene_df.loc[gene_df["products"].apply(has_product), "gene._id"].astype(str))
gene_universe &= protein_genes

print("=== Universe ===")
print(f"Protein-coding proxy gene universe: {len(gene_universe)}")


# ----------------------------
# Curated TF → gene interactions
# ----------------------------
tf_regs = regulon_df.loc[
    regulon_df["regulator._id"].astype(str).str.contains("TFC", na=False),
    ["regulatoryInteractions"],
].copy()

tf_regs = tf_regs[tf_regs["regulatoryInteractions"].apply(is_list)]
tf_regs = tf_regs.explode("regulatoryInteractions", ignore_index=True)

ri_flat = pd.json_normalize(tf_regs["regulatoryInteractions"], sep=".")

ri_gene = ri_flat.loc[
    ri_flat["regulatedEntity.type"].eq("gene"),
    ["regulatedEntity._id"]
].copy()

ri_gene["gene_id"] = ri_gene["regulatedEntity._id"].astype(str)

regulated_genes = set(ri_gene["gene_id"]) & gene_universe

print("\n=== Curated TF regulation ===")
print(f"Genes with ≥1 curated TF interaction: {len(regulated_genes)}")


# ----------------------------
# Build TU → genes map
# ----------------------------
tu_genes = tu_df[["_id", "genes"]].rename(columns={"_id": "tu_id"})
tu_genes = tu_genes.explode("genes", ignore_index=True)

gene_flat = pd.json_normalize(tu_genes["genes"], sep=".")
gene_flat["tu_id"] = tu_genes["tu_id"].astype(str).values
gene_flat = gene_flat.rename(columns={"_id": "gene_id"})

tu_gene_map = gene_flat[["tu_id", "gene_id"]].dropna()
tu_gene_map["gene_id"] = tu_gene_map["gene_id"].astype(str)

# Restrict to protein-coding universe
tu_gene_map = tu_gene_map[tu_gene_map["gene_id"].isin(gene_universe)]


# ----------------------------
# Identify regulated TUs
# ----------------------------
regulated_tus = set(
    tu_gene_map.loc[
        tu_gene_map["gene_id"].isin(regulated_genes),
        "tu_id"
    ]
)

print("\n=== Transcription units ===")
print(f"Total transcription units (with ≥1 protein-coding gene): {tu_gene_map['tu_id'].nunique()}")
print(f"Regulated transcription units: {len(regulated_tus)}")


# ----------------------------
# Genes in regulated TUs
# ----------------------------
genes_in_regulated_tus = set(
    tu_gene_map.loc[
        tu_gene_map["tu_id"].isin(regulated_tus),
        "gene_id"
    ]
)

print("\n=== Gene coverage via regulated TUs ===")
print(f"Genes in ≥1 regulated TU: {len(genes_in_regulated_tus)}")
print(f"Fraction of genes in regulated TUs: {len(genes_in_regulated_tus)/len(gene_universe):.3f}")

genes_not_in_regulated_tus = gene_universe - genes_in_regulated_tus
print(f"Genes NOT in any regulated TU: {len(genes_not_in_regulated_tus)}")
