#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  9 15:00:27 2025

@author: mariepyun
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 12:53:38 2025

@author: brendonmcguinness
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 12 14:09:44 2025
Refactored to isolate COMETS runs via chdir, add error handling,
avoid Decimal keys, load base model once, and clean up workspaces.
"""

import os
import shutil
import tempfile
from copy import deepcopy

import numpy as np
import pandas as pd
import cometspy as c
import cobra
from cobra.io import load_model, read_sbml_model

from helpers_functional import (
    nicheOverlapSaveTimeSeriesMultipleChooseTime,
    fitnessDifferenceSaveTimeSeriesMultipleRandChooseTime,
)

counter = 0

# --- Configuration ---
os.environ['COMETS_HOME'] = '/Applications/comets_macos/comets_2.12.3'
os.environ['GUROBI_COMETS_HOME'] = '/Library/gurobi1202/macos_universal2'
os.environ['GRB_LICENSE_FILE']='/Library/gurobi1202/macos_universal2/gurobi.lic' 

carbon_pairs = [
    #('EX_glc__D_e', 'EX_succ_e'),
    #('EX_glc__D_e', 'EX_fru_e')#,
    ('EX_glc__D_e', 'EX_cit_e')#,
    #('EX_glc__D_e', 'EX_xyl__D_e')#,

]


#Add the network pairs we will shuffle through
network_pairs = [#('network_files/senterica.xml', 'network_files/vcholerae.xml')#,
                 #('network_files/senterica.xml', 'network_files/sarizonae.xml')#,
                 ('network_files/senterica.xml','network_files/ecolik12.xml'),
                #('network_files/senterica.xml','network_files/saureus.xml')
                ]


ko_bounds   = np.arange(-10, 1)  # -10, -9, …, 0
#ko_bounds   = np.array([-10,-3,-2,-1,0])
M           = 2
s1_conc     = np.linspace(0.005, 0.05, M)
s2_conc     = s1_conc[::-1]
#exchange_rxns = ["EX_glc__D_e", "EX_fru_e", "EX_succ_e", "EX_cit_e", "EX_xyl__D_e"]
exchange_rxns = [
    "EX_glc__D_e", "EX_fru_e", "EX_gal_e", "EX_man_e", "EX_xyl__D_e", 
    "EX_arab__L_e", "EX_ac_e", "EX_lac__D_e", "EX_pyr_e", "EX_succ_e", 
    "EX_mal__L_e", "EX_fum_e", "EX_glyc_e", "EX_cit_e", "EX_akg_e",
    "EX_eth_e", "EX_for_e"
]
output_dir = 'coex_data_mut1_vs_mut2_refactored_test_wCross'
os.makedirs(output_dir, exist_ok=True)

# --- Helpers ---

def make_key(network1, network2, s1_name, s2_name, s1, ko1, ko2, counter): #Added "network1" and "network2" to the key, to better identify which experiment is which
    return (
        network1,
        network2,
        (s1_name, s2_name),
        f"{s1:.4f}",
        ko1,
        ko2,
        counter
    )

def save_results_to_csv(results, fname): #Added some columns to the table so that I better know which line corresponds to which experiment
    rows = []
    for key, vals in results.items():
        (model1, model2, pair, conc_str, ko1, ko2, counter) = key
        src1, src2 = pair
        for idx, v in enumerate(vals):
            rows.append({
                "Network 1": model1, #Added network names as first two columns
                "Network 2": model2,
                "Carbon Source 1": src1,
                "Carbon Source 2": src2,
                "Concentration of CS1": float(conc_str), #Specified that we will be referring to carbon source 1 when we write the concentration (concentraton of cs2 can be deduced as they add up to 0.055)
                "KO Bound Source 1 on Network 1": ko1, #Specified that I will be writing the ko bounds on network1 (not network2, but the bounds on network2 can be inferred. For ex: if we have that network 1 has bounds -10, -3, then we know network 2 has bound -3, -10. Just by design of the setup_mutants() function.)
                "KO Bound Source 2 on Network 2": ko2, #assume if i does not equal j (network 1 on carbon 2) then -10
                "KO Strain Index": (idx % 2) + 1,
                "Value": v,
                "Trial #": counter
            })
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, fname), index=False)
    print(f"✅ Saved: {fname}")

def setup_mutants(base1, base2, net1, net2, s1, s2, ko1, ko2, cross_bound, rxns): #Added the possibility of adding 2 different base models (rather than just one)
    # "base1" and "base2" are the actual loaded models
    # "net1" and "net2" are strings, they are the paths to the network files, e.g. 'network_files/senterica.xml'
    mut1 = deepcopy(base1) 
    mut2 = deepcopy(base2)
    for r in rxns:
        mut1.change_bounds(r, cross_bound, 1000)
        mut2.change_bounds(r, cross_bound, 1000)
    mut1_name = net1.split('/')[-1].split('.')[0] #Remove all the extra stuff from the network file path, and only keep the file name
    mut2_name = net2.split('/')[-1].split('.')[0] #Remove all the extra stuff from the network file path, and only keep the file name
    mut1.id = f"{mut1_name}_{s1}_KO_{ko1}"
    mut1.change_bounds(s1, ko1, 1000)
    mut1.change_bounds(s2, -10, 1000)
    mut2.id = f"{mut2_name}_{s2}_KO_{ko2}"
    mut2.change_bounds(s2, ko2, 1000)
    mut2.change_bounds(s1, -10, 1000)
    return mut1, mut2, mut1.id, mut2.id

def run_in_temp(fn, *args, **kwargs):
    """chdir into a temp dir, run `fn`, then clean up."""
    ws = tempfile.mkdtemp(prefix="comets_run_")
    cwd = os.getcwd()
    try:
        os.chdir(ws)
        return fn(*args, **kwargs)
    finally:
        os.chdir(cwd)
        shutil.rmtree(ws, ignore_errors=True)

# --- Load base model once ---
#base_model = c.model(load_model("iJO1366"))

# --- Storage dicts ---
nd, fd, coex, win, mgg, mgl, sc = (
    {}, {}, {}, {}, {}, {}, {}
)


# --- Main loops ---
for network1, network2 in network_pairs:
    
    #load models for each network
    base_model_1 = c.model(read_sbml_model(network1))
    base_model_2 = c.model(read_sbml_model(network2))
    
    for s1_name, s2_name in carbon_pairs:
        
        for ko in ko_bounds:
            
            for i in range(2): #This for loop will run twice, making different modifications to the networks each time
                
                if i == 0: #Set up mutants one way
                    mut1, mut2, mut1_id, mut2_id = setup_mutants(base_model_1, base_model_2, network1, network2, s1_name, s2_name, ko, ko, -10, exchange_rxns)
                
                if i == 1: #Set up mutant the other way
                    mut1, mut2, mut1_id, mut2_id = setup_mutants(base_model_2, base_model_1, network2, network1, s1_name, s2_name, ko, ko, -10, exchange_rxns)
                           
                    
                for conc1, conc2 in zip(s1_conc, s2_conc):
                    
                    counter += 1 
                    
                    if i == 0: #Make the keys correctly (recall that for the key, I am writing ko1 and ko2 bounds based on the uptake capacities of s. enterica)
                        key = make_key(mut1_id, mut2_id, s1_name, s2_name, conc1, ko, -10, counter) #If i == 0, then s. enterica has uptake capacity ko for carbon source 1 and -10 for carbon source 2, so write the key accordingly
                    
                    if i == 1:
                        key = make_key(mut1_id, mut2_id, s1_name, s2_name, conc1, -10, ko, counter) #If i == 0, then s. enterica has uptake capacity ko for carbon source 2 and -10 for carbon source 1, so write the key accordingly
        
        
                    # Niche
                    try:
                        #Always add S. enterica to the test tube first
                        if i == 0:
                            no, co_val, b1, m1, winner, sc_val = run_in_temp(
                                nicheOverlapSaveTimeSeriesMultipleChooseTime,
                                mut1, mut2, s1_name, s2_name, conc1, conc2, max_cyc=320
                            )
                            
                        if i ==1:
                             no, co_val, b1, m1, winner, sc_val = run_in_temp(
                                 nicheOverlapSaveTimeSeriesMultipleChooseTime,
                                 mut2, mut1, s1_name, s2_name, conc1, conc2, max_cyc=320
                             )
                        
                    except Exception as e:
                        print(f"⚠️  Niche failed {key}: {e}")
                        continue
    
    
                    # Fitness
                    try:
                         #Always add S. enterica to the test tube first
                         if i == 0:
                             fd_val, b3, m3, mg_grad, mg_slope = run_in_temp(
                                 fitnessDifferenceSaveTimeSeriesMultipleRandChooseTime,
                                 mut1, mut2, s1_name, s2_name, conc1, conc2, max_cyc=320
                             )
                         
                         if i == 1:
                             fd_val, b3, m3, mg_grad, mg_slope = run_in_temp(
                                 fitnessDifferenceSaveTimeSeriesMultipleRandChooseTime,
                                 mut2, mut1, s1_name, s2_name, conc1, conc2, max_cyc=320
                             )

                    except Exception as e:
                        print(f"⚠️  Fitness failed {key}: {e}")
                        continue
        
                    nd.setdefault(key, []).append(no)
                    coex.setdefault(key, []).append(co_val)
                    win.setdefault(key, []).append(winner)
                    sc.setdefault(key, []).append(sc_val)
                    fd.setdefault(key, []).append(fd_val)
                    mgg.setdefault(key, []).append(mg_grad)
                    mgl.setdefault(key, []).append(mg_slope)
                    #print(key)

# --- Convert & save ---
for dct, fname in [
    (nd,   "niche_differences.csv"),
    (fd,   "fitness_differences.csv"),
    (coex, "coexistence_results.csv"),
    (win,  "winner_results.csv"),
    (mgg,  "mgg_diff_results.csv"),
    (mgl,  "mgl_diff_results.csv"),
    (sc,   "sc_diff_results.csv"),
]:
    
    for k in dct:
        dct[k] = np.array(dct[k])
        
    save_results_to_csv(dct, fname)