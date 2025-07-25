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
        "KO Bound Source 1", "KO Bound Source 2",
        "Concentration"
    ]

    # load & rename
    data = {}
    for key, (fname, colname) in files.items():
        path = os.path.join(output_dir, fname)
        if not os.path.isfile(path):
            raise FileNotFoundError(f"'{fname}' not found in {output_dir}")
        df = pd.read_csv(path)
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
    out['ND_norm'] = out[niche_col] / (1.0 + out[niche_col])
    
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
    
    # 2) Compute your standard transforms
    df_t = transform_differences(df_pair,
                                 niche_col=x_trait,
                                 fitness_col=y_trait,
                                 base=10)
    df_t['Coex_bool'] = df_t['Coexistence Strength'] > 0

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
    ymin, ymax = -0.1, 2.3

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

    # now overplot your boolean points
    for flag, col, lbl in [(True,'red','Coexist'), (False,'blue','Exclusion')]:
        sub = df_t[df_t['Coex_bool'] == flag]
        ax.scatter(sub['ND_norm_stretched'], sub['FD_ratio'],
                   c=col, s=80, edgecolor='k', label=lbl)

    # 6) overplot the transformed boundaries themselves (optional)
    ax.plot(x_band, y_upper, 'k--', lw=1.5)
    ax.plot(x_band, y_lower, 'k--', lw=1.5)

    # formatting
    ax.set(xlabel='Transformed Niche Difference (ND\')',
           ylabel='Transformed Fitness Difference (FD\')',
           title=f"{network1} vs {network2}\n{cs1.replace('EX_','')} vs {cs2.replace('EX_','')}")
    ax.set_xlim(x_band.min(), x_band.max())
    ax.set_ylim(ymin,      ymax)
    ax.legend(loc='upper left')
    ax.grid(linestyle='--', alpha=0.4)
    plt.tight_layout()
    plt.show()
output_dir = ''
df_all = load_and_merge(output_dir)
cs1, cs2 = "EX_glc__D_e", "EX_cit_e"
for net1, net2 in df_all[["Network 1","Network 2"]].drop_duplicates().values:
    plot_network_boolean_coexistence3(df_all, net1, net2, cs1, cs2)
