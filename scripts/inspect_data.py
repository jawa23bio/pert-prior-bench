"""Look inside an .h5ad before running the benchmark, so you can set the flags.

    python scripts/inspect_data.py --h5ad data/norman.h5ad

Prints obs columns with example values, the most common labels (to find the
control label and the separator used for doubles), and var names / columns
(to check whether genes are symbols or Ensembl IDs).
"""
import argparse

import anndata as ad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    a = ap.parse_args()
    adata = ad.read_h5ad(a.h5ad, backed="r")
    print(adata)
    print("\nobs columns:")
    for c in adata.obs.columns:
        vals = adata.obs[c].astype(str)
        print(f"  {c:30s} {vals.nunique():>7} unique  e.g. {list(vals.unique()[:4])}")
    print("\nvar_names:", list(adata.var_names[:8]))
    print("var columns:", list(adata.var.columns))
    for c in adata.obs.columns:
        vals = adata.obs[c].astype(str)
        if 50 < vals.nunique() < 1000:
            print(f"\nmost common labels in '{c}' (look for the control and the double separator):")
            print(vals.value_counts().head(8).to_string())
            break


if __name__ == "__main__":
    main()
