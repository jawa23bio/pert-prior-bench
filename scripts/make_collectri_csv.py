"""Turn the CollecTRI download from OmniPath into the source,target,weight CSV
the benchmark reads (data/collectri.csv).

Download (TSV, about 4 MB):
    https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes
saved as data/collectri.tsv.

Rules:
- AP-1 and NF-kB complexes (e.g. FOS_JUN) are split into their member genes,
  each inheriting the complex's edges.
- weight = -1 if the edge is inhibition-only, otherwise +1.
- exact duplicate rows are dropped.

    python scripts/make_collectri_csv.py --tsv data/collectri.tsv --out data/collectri.csv
"""
import argparse

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--tsv", default="data/collectri.tsv")
    ap.add_argument("--out", default="data/collectri.csv")
    a = ap.parse_args()

    t = pd.read_csv(a.tsv, sep="\t")
    rows = []
    cols = ["source_genesymbol", "target_genesymbol", "is_stimulation", "is_inhibition"]
    for src, tgt, stim, inh in t[cols].itertuples(index=False):
        w = -1 if (inh and not stim) else 1
        for g in str(src).replace("COMPLEX:", "").split("_"):
            rows.append((g, tgt, w))
    net = pd.DataFrame(rows, columns=["source", "target", "weight"]).drop_duplicates()
    net.to_csv(a.out, index=False)
    print(f"wrote {len(net)} edges, {net['source'].nunique()} TFs to {a.out}")


if __name__ == "__main__":
    main()
