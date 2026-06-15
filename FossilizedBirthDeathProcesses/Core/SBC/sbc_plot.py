#!/usr/bin/env python3
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: sbc_plot.py <ranks.tsv> [out.png] [n_bins]")
    path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else path.replace("_ranks.tsv", "") + "_sbc.png"
    nbins = int(sys.argv[3]) if len(sys.argv) > 3 else 20

    names = open(path).readline().rstrip("\n").split("\t")
    data = np.loadtxt(path, skiprows=1, ndmin=2)
    if data.shape[1] != len(names):
        sys.exit("column count mismatch in " + path)
    R = data.shape[0]

    fig, axes = plt.subplots(1, len(names), figsize=(4 * len(names), 3.4), squeeze=False)
    for j, name in enumerate(names):
        r = data[:, j]
        ax = axes[0][j]
        counts, edges = np.histogram(r, bins=nbins, range=(0.0, 1.0))
        centers = 0.5 * (edges[:-1] + edges[1:])
        expected = R / nbins
        lo, hi = stats.binom.ppf([0.025, 0.975], R, 1.0 / nbins)
        ax.bar(centers, counts, width=1.0 / nbins, color="#6699cc", edgecolor="white")
        ax.axhline(expected, color="black", lw=1)
        ax.axhspan(lo, hi, color="grey", alpha=0.25, lw=0)
        D, p = stats.kstest(r, "uniform")
        ax.set_title("%s\nKS D=%.3f p=%.3f (R=%d)" % (name, D, p, R))
        ax.set_xlabel("normalized rank")
        ax.set_xlim(0, 1)
    axes[0][0].set_ylabel("count")
    fig.tight_layout()
    fig.savefig(out, dpi=300)
    print("wrote", out)


if __name__ == "__main__":
    main()
