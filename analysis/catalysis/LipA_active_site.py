'''Analysis functions that assist with using MD to study polymer-enzyme systems,
   Specifically for the active site of LipA. Will be extended to multiple enzyme
   types in the future.'''

__author__ = 'Joe Laforet Jr.'
__email__ = 'jola3134@colorado.edu'


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from tqdm import tqdm
import os
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import matplotlib.pyplot as plt
import numpy as np

###########################################################################################################################################

# Function to get atom selections based on the provided dictionary

# The dictionary is a mapper between a human-readable string code to a 
# MDAnalysis atom selection
# selection_dict = {
#     "Ser:H": "resid 77 and type H04",
#     "Ser:O": "resid 77 and type OG",
#     "His:N": "resid 156 and type NE2", # nitrogen without the H
#     "His:Ndonor": "resid 156 and type ND1",
#     "His:H": "resid 156 and type H08", # H on the nitrogen, interacts with the Ser-OH, called H07 in the incorrect protonation state, H08 in correct state
#     "Asp:OD1": "resid 133 and type OD1",
#     "Asp:OD2": "resid 133 and type OD2"
# }

def get_selections(universe: mda.Universe,
                   selection_dict: dict):

    selections = {}
    for label, sel_str in selection_dict.items():
        selections[label] = universe.select_atoms(sel_str)
    return selections

###########################################################################################################################################

# NOTE on this function: right now it is hard coded for the specific atom types
# of the catalytic triad Ser-His-Asp. This works for LipA, but you should
# verify that this is the correct triad for your enzyme.
# We also compute the distance to the O as the distance to the midpoint of the 
# two O atoms. This assumption can be changed, because due to resonance the H-bond
# can form between either O of the Asp.

from typing import Dict, List, Tuple
import MDAnalysis as mda

def analyze_distances(universe: mda.Universe,
                      selection_dict: Dict[str, str],
                      interaction_pairs: List[Tuple[str, str]],
                      eq_offset: int = 500) -> Dict:
    """
    Analyze distances between specified atom selection pairs in a molecular dynamics simulation.

    Parameters
    ----------
    universe : mda.Universe
        MDAnalysis Universe object containing trajectory and topology data.

    selection_dict : dict[str, str]
        Dictionary mapping selection names (e.g., "resA") to valid MDAnalysis atom selection strings.

    interaction_pairs : list[tuple[str, str]]
        List of tuples specifying which pairs of selection names (keys from `selection_dict`)
        to compute distances between. Each tuple should be of the form ("sel1", "sel2").

    eq_offset : int, default=500
        Number of initial frames to skip when analyzing distances (used to exclude equilibration).

    Returns
    -------
    dict
        Dictionary mapping each selection pair (e.g., "resA-resB") to a 2D numpy array of
        distances for each frame (shape: [n_frames - eq_offset, n_pairs_per_frame]).
    """

    selections = get_selections(universe, selection_dict)
    n_frames = len(universe.trajectory[eq_offset:])
    
    all_distances = {f"{label1}-{label2}": np.zeros(n_frames) for label1, label2 in interaction_pairs}
    all_distances["His:H-Asp:Oavg"] = np.zeros(n_frames)

    # Get indices for all atoms at once for fast access
    atom_indices = {label: sel.indices for label, sel in selections.items()}

    for i, ts in enumerate(universe.trajectory[eq_offset:]):
        for label1, label2 in interaction_pairs:
            if (label1 == "His:H" and label2 in ["Asp:OD1", "Asp:OD2"]) or label2 == "Asp:Oavg":
                continue  # Skip individual His-H to Asp-ODX
            
            pos1 = universe.atoms[atom_indices[label1]].positions
            pos2 = universe.atoms[atom_indices[label2]].positions
            dist = distances.distance_array(pos1, pos2)[-1]  # get the last value (usually a single atom)
            all_distances[f"{label1}-{label2}"][i] = dist
        
        # Vectorized average His:H to Asp:OD1 and Asp:OD2

        pos_his = universe.atoms[atom_indices["His:H"]].positions

        pos_od1 = universe.atoms[atom_indices["Asp:OD1"]].positions
        pos_od2 = universe.atoms[atom_indices["Asp:OD2"]].positions
        
        d1 = distances.distance_array(pos_his, pos_od1)[-1]
        d2 = distances.distance_array(pos_his, pos_od2)[-1]
        
        all_distances["His:H-Asp:Oavg"][i] = 0.5 * (d1 + d2)

    return all_distances

###########################################################################################################################################

# New function to compare interactions at different temperatures on the same plot
def compare_interactions(u_list, temp_reps, selection_dict,
                         plotting_pairs, title, eq_offset=500):
    start_index = 0
    temperature_data = {label: {} for label, _ in plotting_pairs}
    
    for temp, reps in temp_reps.items():
        for label, _ in plotting_pairs:
            split_label = label.split("-")
            label1 = split_label[0]
            label2 = split_label[1]
            combined_dists = []
            
            for i in range(reps):
                u = u_list[start_index + i][0]
                dists = analyze_distances(u, selection_dict, [(label1, label2)], eq_offset=eq_offset)[f"{label1}-{label2}"]
                print(f"Length of distances: {len(dists)}")
                combined_dists.extend(dists.flatten())
            
            temperature_data[label][temp] = np.array(combined_dists)
        start_index += reps
    
    # Plot the histogram comparing interactions at different temperatures
    plt.figure(figsize=(10, 6))
    colors = ['blue', 'orange', 'green', 'red', 'yellow', 'purple']
    
    all_dists = []

    for idx, (label, temp) in enumerate(plotting_pairs):
        dists = temperature_data[label][temp]
        n, bins, _ = plt.hist(dists, bins=25, density=True, alpha=0.75, color=colors[idx % len(colors)], label=f'{label} at {temp}K')
        all_dists.append((f'{label} at {temp}K', n, bins, dists))
    
    plt.legend()
    plt.xlabel('Interatomic Distance (Å)')
    plt.ylabel('Probability')
    plt.title(f'{title}: Interatomic distance of catalytic atoms')
    plt.savefig(f"{title}.svg")
    plt.show()
    return all_dists
