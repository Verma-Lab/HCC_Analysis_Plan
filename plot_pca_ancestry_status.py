"""
PCA plot colored by pre-assigned ancestry label, with case/control status
shown by marker shape.

This script does NOT plot 1KG individuals
as a visual reference overlay. It plots only study participants, using
genetically inferred ancestry labels already present in the covariate file
(e.g., EUR, AFR, AMR, EAS, SAS).

Color = ancestry group
Marker = case (X) vs control (circle)

--- Mode 1: covariate file already has PC columns (no --evec needed) ---

    python plot_pca_ancestry_status.py \
        --covar covariates_hcc_cancer_free_controls_exomepcs_ancestry_column_3.0.csv \
        --id_col person_id \
        --ancestry_col Ancestry \
        --status_col HCC_Cancer_Free \
        --case_value 1 \
        --control_value 0 \
        --pc_x_col PC1 \
        --pc_y_col PC2 \
        --output pmbb_hcc_pca_ancestry.pdf

--- Mode 2: PCs come from a separate .evec file ---

    python plot_pca_ancestry_status.py \
        --evec imputed_PCA_with1KG.evec \
        --covar covariate_file.csv \
        --id_col person_id \
        --ancestry_col Ancestry \
        --status_col HCC_Cancer_Free \
        --case_value 1 \
        --control_value 0 \
        --output pmbb_hcc_pca_ancestry.pdf
"""

import argparse
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# ── Ancestry color settings ────────────────────────────────────────────────────
# Colors match the 1KG superpop style from the HCC paper for consistency
ANCESTRY_COLORS = {
    "AFR": "#B89BC8",   # lavender
    "AMR": "#56A86B",   # green
    "EAS": "#3AADA8",   # teal
    "EUR": "#3F7FBF",   # blue
    "SAS": "#E0BC3C",   # gold
}

# Marker shape by case/control (matches plot_pca_cases_controls.py style)
STATUS_MARKERS = {
    "Case":    ("X",  40, 0.90),   # filled X, red/pink
    "Control": ("+",  30, 0.75),   # plus, orange
}


def parse_args():
    p = argparse.ArgumentParser(
        description="PCA plot colored by ancestry, case vs control by marker shape."
    )
    p.add_argument("--evec", default=None,
                   help="EIGENSOFT .evec file (SampleID PC1 PC2 ... Population). "
                        "Optional if the covariate file already contains PC columns.")
    p.add_argument("--covar", required=True,
                   help="Covariate file (tab or comma separated) with ancestry and "
                        "case/control columns")
    p.add_argument("--id_col", default="IID",
                   help="ID column name in covariate file (default: IID)")
    p.add_argument("--ancestry_col", default="ancestry",
                   help="Ancestry label column in covariate file (default: ancestry). "
                        "Expected values: EUR, AFR, AMR, EAS, SAS (or a subset).")
    p.add_argument("--status_col", default="STATUS",
                   help="Case/control status column (default: STATUS)")
    p.add_argument("--case_value", default="Case",
                   help="Value in status_col that indicates a case (default: Case). "
                        "Can be a string or integer (e.g., 1).")
    p.add_argument("--control_value", default="Control",
                   help="Value in status_col that indicates a control (default: Control). "
                        "Can be a string or integer (e.g., 0).")
    p.add_argument("--sep", default=None,
                   help="Delimiter for covariate file. Auto-detected from extension if "
                        "not specified (tab for .txt, comma for .csv).")
    p.add_argument("--output", default="pca_ancestry_status.pdf",
                   help="Output file path (.pdf / .png / .svg)")
    p.add_argument("--title", default=None,
                   help="Plot title (optional)")
    p.add_argument("--pc_x", type=int, default=1,
                   help="Which PC to use on X axis when reading from .evec (default: 1)")
    p.add_argument("--pc_y", type=int, default=2,
                   help="Which PC to use on Y axis when reading from .evec (default: 2)")
    p.add_argument("--pc_x_col", default=None,
                   help="Column name for X-axis PC in covariate file (e.g. PC1). "
                        "Required when --evec is not provided.")
    p.add_argument("--pc_y_col", default=None,
                   help="Column name for Y-axis PC in covariate file (e.g. PC2). "
                        "Required when --evec is not provided.")
    p.add_argument("--flip_pc1", action="store_true",
                   help="Negate PC1 values")
    p.add_argument("--flip_pc2", action="store_true",
                   help="Negate PC2 values")
    p.add_argument("--controls_first", action="store_true", default=True,
                   help="Plot controls before cases so cases appear on top (default: True)")
    return p.parse_args()


def load_evec(path, pc_x=1, pc_y=2):
    """Parse EIGENSOFT .evec file, return DataFrame with IID, PC_X, PC_Y."""
    rows = []
    with open(path) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if not line:
                continue
            if i == 0:      # eigenvalue header
                continue
            parts = line.split()
            if len(parts) < 3:
                continue
            sample_id = parts[0]
            pc_vals   = parts[1:-1]
            try:
                x = float(pc_vals[pc_x - 1])
                y = float(pc_vals[pc_y - 1])
            except (IndexError, ValueError):
                continue
            rows.append({"IID": sample_id, "PC_X": x, "PC_Y": y})
    df = pd.DataFrame(rows)
    print(f"Loaded {len(df)} individuals from .evec")
    return df


def load_covar(path, sep=None):
    """Load covariate file, auto-detecting delimiter from extension."""
    if sep is None:
        sep = "\t" if path.endswith(".txt") else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = df.columns.str.strip()
    print(f"Loaded covariate file: {len(df)} rows, columns: {list(df.columns)}")
    return df


def main():
    args = parse_args()

    # ── 1. Load covariate file ────────────────────────────────────────────────
    covar = load_covar(args.covar, args.sep)

    for col in [args.id_col, args.ancestry_col, args.status_col]:
        if col not in covar.columns:
            raise ValueError(
                f"Column '{col}' not found in covariate file. "
                f"Available columns: {list(covar.columns)}"
            )

    covar[args.id_col] = covar[args.id_col].astype(str)

    # ── 2. Get PC values — from covariate file or .evec ───────────────────────
    if args.evec is None:
        # Mode 1: PCs are columns in the covariate file
        if args.pc_x_col is None or args.pc_y_col is None:
            raise ValueError(
                "When --evec is not provided, you must specify --pc_x_col and --pc_y_col "
                "pointing to PC columns in the covariate file."
            )
        for col in [args.pc_x_col, args.pc_y_col]:
            if col not in covar.columns:
                raise ValueError(
                    f"PC column '{col}' not found in covariate file. "
                    f"Available columns: {list(covar.columns)}"
                )
        covar = covar.rename(columns={args.pc_x_col: "PC_X", args.pc_y_col: "PC_Y"})
        covar = covar.rename(columns={args.id_col: "IID"})
        df_pcs = covar[["IID", "PC_X", "PC_Y", args.ancestry_col, args.status_col]].copy()
        df_pcs = df_pcs.rename(columns={args.ancestry_col: "ancestry", args.status_col: "_raw_status"})
        if args.flip_pc1:
            df_pcs["PC_X"] = -df_pcs["PC_X"]
        if args.flip_pc2:
            df_pcs["PC_Y"] = -df_pcs["PC_Y"]
        print(f"Using PC columns '{args.pc_x_col}' and '{args.pc_y_col}' from covariate file")
    else:
        # Mode 2: PCs from .evec file, joined to covariate file
        evec = load_evec(args.evec, args.pc_x, args.pc_y)
        if args.flip_pc1:
            evec["PC_X"] = -evec["PC_X"]
        if args.flip_pc2:
            evec["PC_Y"] = -evec["PC_Y"]
        df_pcs = evec.merge(
            covar[[args.id_col, args.ancestry_col, args.status_col]].rename(
                columns={args.id_col: "IID", args.ancestry_col: "ancestry",
                         args.status_col: "_raw_status"}
            ),
            on="IID", how="inner"
        )
        print(f"Matched {len(df_pcs)} individuals (evec ∩ covariate file)")

    # Normalize status values: map case_value → "Case", control_value → "Control"
    def parse_val(v):
        """Try int conversion so --case_value 1 matches integer column."""
        try:
            return int(v)
        except (ValueError, TypeError):
            return v

    case_val    = parse_val(args.case_value)
    control_val = parse_val(args.control_value)

    df_pcs["_status"] = df_pcs["_raw_status"].map(
        {case_val: "Case", control_val: "Control"}
    )
    n_unmapped = df_pcs["_status"].isna().sum()
    if n_unmapped > 0:
        unique_vals = df_pcs["_raw_status"].unique()
        print(f"WARNING: {n_unmapped} rows had status values not matching "
              f"'{case_val}' or '{control_val}'. Unique values: {unique_vals}")

    df = df_pcs.drop(columns=["_raw_status"])

    # Drop rows with missing ancestry or status, and remove UNKNOWN* labels
    before = len(df)
    df = df.dropna(subset=["ancestry", "_status"])
    df = df[~df["ancestry"].str.upper().str.startswith("UNKNOWN")]
    if len(df) < before:
        print(f"  Dropped {before - len(df)} rows with missing/unknown ancestry or status")

    # Summary
    print(f"\nCase/Control breakdown:")
    print(df["_status"].value_counts().to_string())
    print(f"\nAncestry breakdown:")
    print(df["ancestry"].value_counts().to_string())

    # ── 4. Plot ───────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(7, 6))

    ancestries = sorted(df["ancestry"].unique())

    # Plot controls first (behind), then cases (on top)
    plot_order = ["Control", "Case"]

    for status in plot_order:
        marker, size, alpha = STATUS_MARKERS[status]
        zorder = 2 if status == "Control" else 3
        subset_status = df[df["_status"] == status]

        for anc in ancestries:
            subset = subset_status[subset_status["ancestry"] == anc]
            if subset.empty:
                continue
            color = ANCESTRY_COLORS.get(anc, "#888888")   # grey fallback
            if marker == "+":
                # + is a line-only marker: edgecolors sets its visible color
                ax.scatter(
                    subset["PC_X"], subset["PC_Y"],
                    c=color, s=size, alpha=alpha,
                    marker=marker, zorder=zorder,
                    edgecolors=color, linewidths=1.2, rasterized=True
                )
            else:
                # filled markers (X): white border makes them pop
                ax.scatter(
                    subset["PC_X"], subset["PC_Y"],
                    c=color, s=size, alpha=alpha,
                    marker=marker, zorder=zorder,
                    edgecolors="white", linewidths=0.6, rasterized=True
                )

    # ── 5. Legend ─────────────────────────────────────────────────────────────
    handles = []

    # Ancestry entries (colored circle for each)
    for anc in ancestries:
        color = ANCESTRY_COLORS.get(anc, "#888888")
        h = mlines.Line2D([], [], color=color, marker="o", linestyle="None",
                          markersize=6, markerfacecolor=color, label=anc)
        handles.append(h)

    # Separator (blank) between ancestry and status entries
    handles.append(mlines.Line2D([], [], linestyle="None", label=""))

    # Case/Control entries (black markers showing shape)
    for status, (marker, size, alpha) in STATUS_MARKERS.items():
        h = mlines.Line2D([], [], color="black", marker=marker, linestyle="None",
                          markersize=7 if status == "Case" else 6,
                          markerfacecolor="black", label=status + "s")
        handles.append(h)

    ax.legend(handles=handles, title="Ancestry / Status", frameon=True,
              fontsize=8, title_fontsize=9, loc="best")

    # ── 6. Labels & formatting ────────────────────────────────────────────────
    ax.set_xlabel(f"PC{args.pc_x}", fontsize=11)
    ax.set_ylabel(f"PC{args.pc_y}", fontsize=11)
    if args.title:
        ax.set_title(args.title, fontsize=12)

    ax.grid(True, linestyle=":", linewidth=0.4, alpha=0.6)
    plt.tight_layout()
    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    print(f"\nSaved: {args.output}")


if __name__ == "__main__":
    main()

