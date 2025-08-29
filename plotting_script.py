#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 10:57:48 2025

@author: brendonmcguinness
"""
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

def load_and_merge(output_dir):
    """
    Reads the six COMETS CSVs in `output_dir`:
      - niche_differences.csv
      - fitness_differences.csv
      - coexistence_results.csv
      - mgg_diff_results.csv
      - mgl_diff_results.csv
      - sc_diff_results.csv

    Expects each to contain columns:
      ['Network 1','Network 2',
       'Carbon Source 1','Carbon Source 2',
       'KO Bound Source 1','KO Bound Source 2',
       'Concentration','Value', ...]

    Returns a single DataFrame keyed by those seven columns,
    with 'Value' renamed and pivoted into:
      ['Niche Differences','Fitness Differences',
       'Coexistence Strength','MGG Differences',
       'MGL Differences','SC Differences']
    plus a helper column 'Concentration_Carbon_Source_2'.
    """
    # file → desired column name
    files = {
        "niche":   ("niche_differences.csv",    "Niche Differences"),
        "fitness": ("fitness_differences.csv",  "Fitness Differences"),
        "coex":    ("coexistence_results.csv",  "Coexistence Strength"),
        "mgg":     ("mgg_diff_results.csv",     "MGG Differences"),
        "mgl":     ("mgl_diff_results.csv",     "MGL Differences"),
        "sc":      ("sc_diff_results.csv",      "SC Differences")
    }

    merge_keys = [
        "Network 1", "Network 2",
        "Carbon Source 1", "Carbon Source 2",
        "KO Bound Source 1 on Network 1", "KO Bound Source 2 on Network 1",
        "Concentration"
    ]

    # load & rename
    data = {}
    for key, (fname, colname) in files.items():
        path = os.path.join(output_dir, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"'{fname}' not found in {output_dir}")
        df = pd.read_csv(path)
        # 1) concentration
        if "Concentration of CS1" in df.columns:
            df = df.rename(columns={"Concentration of CS1":"Concentration"})
        
        # 2) KO-bound columns
        #    find whatever they were called and force them to your canonical names
        for i in (1, 2):
            pattern = f"KO Bound Source {i} on "
            # look for a column that starts with that
            matches = [c for c in df.columns if c.startswith(pattern)]
            if matches:
                df = df.rename(columns={ matches[0]: f"KO Bound Source {i} on Network 1" })

        df = df.drop(columns=["KO Strain Index"], errors="ignore")
        df["Value"] = pd.to_numeric(df["Value"], errors="coerce")
        # absolute error for niche + fitness
        if key in ("niche", "fitness"):
            df["Value"] = df["Value"].abs()
        df = df.rename(columns={"Value": colname})
        data[key] = df

    # merge all six on the same keys
    df_all = data["niche"]
    for key in ("fitness", "coex", "mgg", "mgl", "sc"):
        df_all = df_all.merge(data[key], on=merge_keys, how="inner")

    # add the mirrored concentration for CS2
    concs = sorted(df_all["Concentration"].unique())
    mirror_map = dict(zip(concs, concs[::-1]))
    df_all["Concentration_Carbon_Source_2"] = df_all["Concentration"].map(mirror_map)

    return df_all



def transform_differences(df, niche_col="Niche Differences", fitness_col="SC Differences", base=10):
    out = df.copy()
    out['FD_ratio'] = np.exp(out[fitness_col])
#    out['FD_ratio'] = 10 ** (out[fitness_col])
    out['ND_raw'] = base ** out[niche_col]
    out['ND_norm'] = np.abs(out[niche_col]) / (1.0 + np.abs(out[niche_col]))
    
    out['ND_norm_stretched'] = out['ND_norm'] / out['ND_norm'].max()
    out['ND_norm_stretched'] = (
         out['ND_norm'] -  out['ND_norm'].min()
    ) / (
         out['ND_norm'].max() -  out['ND_norm'].min()
    )
    return out


def plot_network_boolean_coexistence3(df_all,
                                      network1, network2,
                                      cs1="EX_glc__D_e", cs2="EX_cit_e",
                                      x_trait="Niche Differences",
                                      y_trait="SC Differences"):
    # 1) Subset to this network + carbon pair
    df_pair = df_all.query(
        "`Network 1` == @network1 & `Network 2` == @network2"
        " & `Carbon Source 1` == @cs1 & `Carbon Source 2` == @cs2"
    )
    if df_pair.empty:
        raise ValueError(f"No data for {network1}/{network2} on {cs1} vs {cs2}")
    df_pair = df_all.query(
        "`Network 1` == @net1 & `Network 2` == @net2"
        " & `Carbon Source 1` == @cs1 & `Carbon Source 2` == @cs2"
    )
    print("subset size:", df_all.shape)
    #print(df_all[["Niche Differences","SC Differences"]].describe())

    # 2) Compute your standard transforms
    #df_t = df_all
    df_t = transform_differences(df_all,
                                 niche_col=x_trait,
                                 fitness_col=y_trait,
                                 base=10)
    df_t['Coex_bool'] = df_t['Coexistence Strength'] > 0
    print(df_t[["Niche Differences","SC Differences"]].describe())

    # 3) Build the transformed ±ND boundaries once
    t = np.linspace(0, df_pair[x_trait].max(), 500)
    up   = pd.DataFrame({x_trait: t, y_trait:  t})
    lo   = pd.DataFrame({x_trait: t, y_trait: -t})
    up_t = transform_differences(up,   niche_col=x_trait, fitness_col=y_trait, base=10)
    lo_t = transform_differences(lo,   niche_col=x_trait, fitness_col=y_trait, base=10)
    x_band   = up_t['ND_norm_stretched']
    y_upper  = up_t['FD_ratio']
    y_lower  = lo_t['FD_ratio']

    # 4) Pick your y‐limits (you can also derive these from df_t['FD_ratio'].min()/max())
    ymin, ymax = -0.1, 15

    # 5) Plot!
    fig, ax = plt.subplots(figsize=(7,5))
    # fill the coexistence band in red
    ax.fill_between(x_band, y_lower, y_upper,
                    facecolor='r', alpha=0.2,
                    label='Transformed coexistence region')
    # fill exclusion above it in blue
    ax.fill_between(x_band, y_upper, ymax,
                    facecolor='b', alpha=0.2)
    # fill exclusion below it in blue
    ax.fill_between(x_band, ymin, y_lower,
                    facecolor='b', alpha=0.2)
    print("» total points in df_t:", len(df_t))
    print(df_t[["ND_norm_stretched","FD_ratio","Coex_bool"]])
    # now overplot your boolean points
    for flag, col, lbl in [(True,'red','Coexist'), (False,'blue','Exclusion')]:
        sub = df_t[df_t['Coex_bool'] == flag]
        ax.scatter(sub['ND_norm_stretched'], sub['FD_ratio'],
                   c=col, s=80, edgecolor='k', label=lbl)
        print(f"   Coex_bool={flag:5} → {len(sub):3d} points")
        # show the first few so you can eyeball valid ND/FD
        print(sub[["ND_norm_stretched","FD_ratio"]])
    # 6) overplot the transformed boundaries themselves (optional)
    ax.plot(x_band, y_upper, 'k--', lw=1.5)
    ax.plot(x_band, y_lower, 'k--', lw=1.5)

    # formatting
    ax.set(xlabel='Transformed Niche Difference (ND\')',
           ylabel='Transformed Fitness Difference (FD\')')
    #ax.set_xlim(x_band.min(), x_band.max())
    ax.set_ylim(ymin,ymax)
    # after you create df_t…
    x_min = df_t["ND_norm_stretched"].min()
    x_max = df_t["ND_norm_stretched"].max()
    ax.set_xlim(x_min, x_max)
    ax.relim()
    ax.autoscale_view()
    ax.legend(loc='upper left')
    ax.grid(linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.show()
    
def plot_network_boolean_coexistence_untransformed(
    df_all,
    network1, network2,
    cs1="EX_glc__D_e", cs2="EX_cit_e",
    x_trait="Niche Differences",
    y_trait="SC Differences",
    x_from_zero=True,
    ylim=(-5, 5),
    fit_decision_boundary=True,
    C=1e6,                    # large C ≈ weak regularization
    max_iter=2000,
    savepath=None,
    show=True
):
    """
    Plot 'untransformed' ND vs FD for a given network/carbon-source pair.

    - X axis: absolute(Niche Differences)
    - Y axis: Fitness Differences (sign preserved, e.g., 'SC Differences')
    - Shades coexistence wedge between y = ± x
    - Colors points by Coex_bool (red = coexist, blue = exclusion)
    - Optionally fits and overlays a logistic regression decision boundary

    Parameters
    ----------
    df_all : pd.DataFrame
        Output of `load_and_merge(...)` with columns including:
        ['Niche Differences','SC Differences','Coexistence Strength',
         'Network 1','Network 2','Carbon Source 1','Carbon Source 2', ...]
    network1, network2 : str
        Network identifiers to subset.
    cs1, cs2 : str
        Carbon sources to subset (order matters: “CS1 vs CS2”).
    x_trait : str
        Column name for niche differences (default: 'Niche Differences').
    y_trait : str
        Column name for fitness differences (default: 'SC Differences').
    x_from_zero : bool
        If True, draw wedge from x=0 to max(|ND|); else from data min..max.
    ylim : (float, float)
        Y-axis limits.
    fit_decision_boundary : bool
        If True, fit LogisticRegression and overlay boundary.
    C, max_iter : numeric, int
        LogisticRegression hyperparameters.
    savepath : str or None
        If provided, saves the figure to this path.
    show : bool
        If True, calls plt.show().

    Returns
    -------
    result : dict
        Keys: 'clf', 'slope', 'intercept', 'xlim', 'ylim', 'subset'
        (clf/slope/intercept are None if boundary not fit or only one class)
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # Subset to the requested pair
    q = (
        (df_all["Network 1"] == network1) &
        (df_all["Network 2"] == network2) &
        (df_all["Carbon Source 1"] == cs1) &
        (df_all["Carbon Source 2"] == cs2)
    )
    df = df_all.loc[q].copy()
    if df.empty:
        raise ValueError(f"No data for ({network1}, {network2}) on {cs1} vs {cs2}")

    # Labels / features
    df["ND_abs"] = df[x_trait].abs()
    df["FD"] = pd.to_numeric(df[y_trait], errors="coerce")
    df["Coex_bool"] = df["Coexistence Strength"] > 0.0
    df = df.dropna(subset=["ND_abs", "FD", "Coexistence Strength"])

    # X range for wedge/boundaries
    if x_from_zero:
        x_min = 0.0
        x_max = float(df["ND_abs"].max()) if df["ND_abs"].size else 1.0
    else:
        x_min = float(df["ND_abs"].min())
        x_max = float(df["ND_abs"].max())

    x_vals = np.linspace(x_min, x_max, 500)
    y_upper = +x_vals
    y_lower = -x_vals

    # Prepare figure
    fig, ax = plt.subplots(figsize=(7, 5))

    # Fill regions
    y_min, y_max = ylim
    ax.fill_between(x_vals, y_lower, y_upper, facecolor='r', alpha=0.2, label='Coexistence')
    ax.fill_between(x_vals, y_upper, y_max, facecolor='b', alpha=0.2)
    ax.fill_between(x_vals, y_min, y_lower, facecolor='b', alpha=0.2)

    # Draw FD = ± ND lines
    ax.plot(x_vals, y_upper, 'k-', lw=2, label='FD = ND')
    ax.plot(x_vals, y_lower, 'k-', lw=2)

    # Scatter points colored by coexistence boolean
    for flag, col, lbl in [(True, 'red', 'Coexist'), (False, 'blue', 'Exclusion')]:
        sub = df.loc[df["Coex_bool"] == flag]
        ax.scatter(sub["ND_abs"], sub["FD"], c=col, s=80, edgecolor='k', label=lbl, alpha=0.9)

    # Optional: logistic boundary
    slope = intercept = None
    clf = None
    if fit_decision_boundary and df["Coex_bool"].nunique() == 2:
        try:
            from sklearn.linear_model import LogisticRegression
            X = df[["ND_abs", "FD"]].to_numpy()
            y = df["Coex_bool"].astype(int).to_numpy()
            clf = LogisticRegression(
                penalty="l2",
                C=C,
                solver="lbfgs",
                tol=1e-4,
                max_iter=max_iter
            ).fit(X, y)

            w0, w1 = clf.coef_[0]
            b = clf.intercept_[0]

            # If w1 ≈ 0, the boundary is nearly vertical: w0*x + b = 0
            eps = 1e-12
            if abs(w1) < eps:
                # vertical line at x = -b / w0
                if abs(w0) > eps:
                    x_v = -b / w0
                    ax.axvline(x_v, linestyle='--', color='k', lw=2, label='Decision boundary')
                # else degenerate; skip drawing
            else:
                # y = -(w0*x + b)/w1
                y_dec = -(w0 * x_vals + b) / w1
                #ax.plot(x_vals, y_dec, 'k--', lw=2, label='Decision boundary')
                slope = -w0 / w1
                intercept = -b / w1
                # Optional: show equation text
                # ax.text(x_min, y_max*0.9, f"y = {slope:.3f} x + {intercept:.3f}")

        except Exception as e:
            print(f"[warn] Logistic boundary not drawn: {e}")

    # Axes / formatting
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_xlabel('Raw Niche Difference, ND', fontsize=13)
    ax.set_ylabel('Raw Fitness Difference, FD', fontsize=13)
    ax.legend(loc='upper left', fontsize=11)
    ax.grid(linestyle='--', alpha=0.35)
    plt.tight_layout()

    if savepath:
        plt.savefig(savepath, bbox_inches='tight')

    if show:
        plt.show()
    else:
        plt.close(fig)

    return None 
    """{
        "clf": clf,
        "slope": slope,
        "intercept": intercept,
        "xlim": (x_min, x_max),
        "ylim": (y_min, y_max),
        "subset": df.copy()
    }"""


def plot_ko1_vs_niche_for(
    df_all,
    network1, network2,
    cs1, cs2,
    concentration,
    color_by=None,     # None | 'Coex_bool' | 'Coexistence Strength' (or any numeric column)
    s=80
):
    """
    Scatter: (KO Bound Source 1 on Network 1) vs (Niche Differences),
    filtered to a specific network pair, carbon sources, and concentration.
    """
    # Ensure numeric concentration for robust matching
    df = df_all.copy()
    df["Concentration"] = pd.to_numeric(df["Concentration"], errors="coerce")

    mask = (
        (df["Network 1"] == network1) &
        (df["Network 2"] == network2) &
        (df["Carbon Source 1"] == cs1) &
        (df["Carbon Source 2"] == cs2) &
        (np.isfinite(df["Concentration"])) &
        (np.isclose(df["Concentration"], float(concentration), rtol=0.0, atol=1e-12))
    )
    sub = df.loc[mask, [
        "KO Bound Source 1 on Network 1",
        "Niche Differences",
        "Coexistence Strength"
    ]].copy()

    if sub.empty:
        raise ValueError(
            f"No rows for {network1}/{network2} on {cs1} vs {cs2} at concentration={concentration}."
        )

    # Optional coloring choices
    c_kw = {}
    if color_by == "Coex_bool":
        sub["Coex_bool"] = sub["Coexistence Strength"] > 0
        colors = sub["Coex_bool"].map({True: "red", False: "blue"}).values
        c_kw = {"c": colors}
    elif color_by and color_by in df.columns:
        # color by a numeric column (e.g., 'Coexistence Strength')
        vals = df.loc[mask, color_by].to_numpy()
        c_kw = {"c": vals, "cmap": "viridis"}
    # else: default matplotlib color

    plt.figure(figsize=(7,5))
    plt.scatter(
        sub["KO Bound Source 1 on Network 1"],
        sub["Niche Differences"],
        edgecolor="k",
        alpha=0.9,
        s=s,
        **c_kw
    )
    plt.xlabel("KO Bound Source 1 on Network 1", fontsize=12)
    plt.ylabel("Niche Differences", fontsize=12)
    plt.title(
        f"{network1} vs {network2} | {cs1} vs {cs2} | Conc={concentration}",
        fontsize=12
    )
    if color_by == "Coex_bool":
        from matplotlib.lines import Line2D
        legend_elems = [
            Line2D([0],[0], marker='o', color='w', markerfacecolor='red', markeredgecolor='k', label='Coexist', markersize=8),
            Line2D([0],[0], marker='o', color='w', markerfacecolor='blue', markeredgecolor='k', label='Exclusion', markersize=8),
        ]
        plt.legend(handles=legend_elems, loc='best')
    elif color_by:
        plt.colorbar(label=color_by)
    plt.tight_layout()
    plt.show()


output_dir = 'cobra_results' #'coex_data_marie'
df_all = load_and_merge(output_dir)
cs1, cs2 = "EX_glc__D_e", "EX_cit_e"
net1 = "salmonella" #"senterica" #"_EX_glc__D_e_KO_-10"
net2 = "iJO1366" #"ecolik12" #"_EX_cit_e_KO_-10"
plot_network_boolean_coexistence3(df_all, net1, net2, cs1, cs2)
#for net1, net2 in df_all[["Network 1","Network 2"]].drop_duplicates().values:
#    plot_network_boolean_coexistence3(df_all, net1, net2, cs1, cs2)
untrans_bound = 8
_ = plot_network_boolean_coexistence_untransformed(
        df_all, net1, net2, cs1, cs2,
        x_trait="Niche Differences",
        y_trait="SC Differences",
        ylim=(-untrans_bound, untrans_bound),
        fit_decision_boundary=True,
        savepath="untransformed_glc_cit.pdf"
    )

plot_ko1_vs_niche_for(df_all, net1, net2, cs1, cs2, 0.05)