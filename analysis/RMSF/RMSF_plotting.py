'''Functions used to plot the RMSF of various portions of the protein in a simulation'''

__author__ = 'Joe Laforet Jr.'
__email__ = 'jola3134@colorado.edu'

import MDAnalysis as mda
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MDAnalysis.analysis import rms, align

from biopandas.pdb import PandasPdb

import json
from collections import defaultdict
from tqdm import tqdm
from matplotlib import cm
from matplotlib.cm import get_cmap

# Suppress some MDAnalysis warnings about writing PDB files
import warnings
warnings.filterwarnings('ignore')

from typing import Union

# Function to read in a JSON file
def read_json_file(file_name):
    with open(file_name, 'r') as f:
        data = json.load(f)
    return data


#########################################################################################

# Example parameters:
# json_files = ['1000ns_noPoly_allResids_rmsf_results.json',          '1000ns_EGMA_allResids_rmsf_results.json',
#              '1000ns_SBMA-EGPMA_100:0_allResids_rmsf_results.json', '1000ns_SBMA-EGPMA_99:1_allResids_rmsf_results.json',
#              '1000ns_SBMA-EGPMA_98:2_allResids_rmsf_results.json',  '1000ns_SBMA-EGPMA_95:5_allResids_rmsf_results.json',
#              '1000ns_SBMA-EGPMA_90:10_allResids_rmsf_results.json']
              
# system_names = ['noPoly', 'EGMA', '100:0', '99:1', '98:2', '95:5', '90:10']
# all_res_ids = np.arange(1,182)

# Used for plotting the RMSF for the whole protein
def plot_per_residue_rmsf_from_json(json_files: list[str],
                                    system_names: list[str],
                                    all_res_ids: list[int],
                                    temperatures: list[str]=["293"], 
                                    ylim: Union[None, int, float]=None, 
                                    highlight_resids: list[str]=None,
                                    second_highlight_resids: list[str]=None,
                                    legend_title: str="SBMA:EGPMA",
                                    out_title: str="Per_Residue_RMSF.pdf"):
    """
    Plots per-residue RMSF (mean of time-resolved RMSD) for one or more temperatures.
    Each temperature gets its own subplot (row). Contiguous highlight regions are shaded.

    Parameters:
        - json_files: list of file paths
        - system_names: list of system labels
        - all_res_ids: array of all residue IDs to use (same across all files)
        - temperatures: str or list of temperature keys (e.g., ["293", "310"])
        - highlight_resids: list of residue indices (int or str) to highlight as shaded regions

    Returns:
        results[temp][system_name] = {
            'resids': [...],
            'per_rep_rmsfs': ndarray[num_reps, num_residues]
        }
    """
    if isinstance(temperatures, str) or isinstance(temperatures, int):
        temperatures = [str(temperatures)]
    else:
        temperatures = [str(t) for t in temperatures]

    tab10 = get_cmap("tab10")
    system_colors = [tab10(i % 10) for i in range(len(system_names))]
    all_per_rep_rmsfs = defaultdict(dict)

    legend_entries = []

    # Standardize highlight residue list to sorted integers
    highlight_resids = sorted(set(int(r) for r in highlight_resids)) if highlight_resids else []

    second_highlight_resids = sorted(set(int(r) for r in second_highlight_resids)) if second_highlight_resids else []

    # Group consecutive residues into blocks
    def group_consecutive(resid_list):
        groups = []
        if not resid_list:
            return groups
        start = prev = resid_list[0]
        for r in resid_list[1:]:
            if r == prev + 1:
                prev = r
            else:
                groups.append((start, prev))
                start = prev = r
        groups.append((start, prev))
        return groups

    highlight_blocks = group_consecutive(highlight_resids)
    second_highlight_blocks = group_consecutive(second_highlight_resids)

    n_temps = len(temperatures)
    fig, axes = plt.subplots(n_temps, 1, figsize=(14, 4 * n_temps), sharex=True)
    if n_temps == 1:
        axes = [axes]  # Make iterable

    for ax, temp in zip(axes, temperatures):
        for file_path, system_name, color in zip(json_files, system_names, system_colors):
            data_dict = read_json_file(file_path)

            if temp not in data_dict:
                print(f"Temperature {temp} not found in {file_path}, skipping.")
                continue

            temp_data = data_dict[temp]
            resids = all_res_ids.astype(str)
            num_residues = len(resids)

            per_rep_rmsfs = []

            # Get number of replicates from first residue
            first_res = resids[0]
            print(temp_data.keys())
            print(file_path)
            num_reps = len(temp_data[first_res])

            for rep_idx in range(num_reps):
                rep_rmsf = []
                for resid in resids:
                    rmsd_vals = temp_data[resid][rep_idx]['rmsd']
                    rep_rmsf.append(np.mean(rmsd_vals))
                per_rep_rmsfs.append(rep_rmsf)

            per_rep_rmsfs = np.array(per_rep_rmsfs)
            avg_rmsf = np.mean(per_rep_rmsfs, axis=0)
            sem_rmsf = np.std(per_rep_rmsfs, axis=0, ddof=1) / np.sqrt(num_reps)

            all_per_rep_rmsfs[temp][system_name] = {
                'resids': resids,
                'per_rep_rmsfs': per_rep_rmsfs
            }

            resids_int = [int(r) for r in resids]
            line, = ax.plot(resids_int, avg_rmsf, label=system_name, color=color, linewidth=2)
            ax.fill_between(resids_int, avg_rmsf - sem_rmsf, avg_rmsf + sem_rmsf, color=color, alpha=0.3)

            # Only collect legend entries once (first temp loop)
            if temp == temperatures[0]:
                legend_entries.append((line, system_name))

        # Step 1: Track yellow-shaded residue indices
        yellow_indices = set()
        for start, end in highlight_blocks:
            ax.axvspan(start - 0.5, end + 0.5, color='yellow', alpha=0.2)
            yellow_indices.update(range(start, end + 1))

        # Step 2: Shade blue only for non-overlapping regions
        for start, end in second_highlight_blocks:
            # Compute set of this blue block's indices
            blue_indices = set(range(start, end + 1))
            non_overlap = blue_indices - yellow_indices

            if non_overlap:
                # Find contiguous blocks within non-overlapping region
                sorted_indices = sorted(non_overlap)
                sub_start = sorted_indices[0]
                for i in range(1, len(sorted_indices)):
                    if sorted_indices[i] != sorted_indices[i - 1] + 1:
                        ax.axvspan(sub_start - 0.5, sorted_indices[i - 1] + 0.5, color='cyan', alpha=0.2)
                        sub_start = sorted_indices[i]
                ax.axvspan(sub_start - 0.5, sorted_indices[-1] + 0.5, color='cyan', alpha=0.2)

        if ylim:
            ax.set_ylim(ylim)

        #ax.set_ylabel("Time-averaged RMSD (Å)")
        #ax.set_title(f"Per-Residue RMSF at {int(temp)-273}°C")
        #ax.grid(True, alpha=0.3)
        #ax.legend(title=f"{legend_title}")

    axes[-1].set_xlabel("Residue Number")
    # Shared legend
    handles, labels = zip(*legend_entries)
    fig.legend(handles, labels, title=legend_title, loc='center left', bbox_to_anchor=(1.01, 0.5), frameon=False)
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave room on right for legend
    plt.savefig(out_title)
    plt.show()

    return all_per_rep_rmsfs

##########################################################################################

