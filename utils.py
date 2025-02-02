import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors, PandasTools
from utils import *



def read_smi(path, delimiter='\t', titleLine=False):
    """Read SMILES from .txt, .sdf, or .smi files."""
    result = []
    if path.endswith('.txt'):
        with open(path, 'r') as f:
            for smi in tqdm(f.readlines(), desc='Reading SMILES'):
                if Chem.MolFromSmiles(smi):
                    result.append(smi.strip())
    elif path.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(path)
        for mol in tqdm(supplier, desc='Reading SDF'):
            if mol:
                result.append(Chem.MolToSmiles(mol))
    elif path.endswith('.smi'):
        supplier = Chem.SmilesMolSupplier(path, delimiter=delimiter, titleLine=titleLine)
        for mol in tqdm(supplier, desc='Reading SMILES'):
            if mol:
                result.append(Chem.MolToSmiles(mol))
    return result

def update_dataframe(x, bound):
    """Update dataframe with molecular descriptors based on provided thresholds."""
    if isinstance(x, list):
        x = pd.DataFrame(x, columns=['SMILES'])
        x['ROMol'] = x['SMILES'].apply(Chem.MolFromSmiles)
    for key, _ in bound.items():
        for desc, func in Descriptors._descList:
            if key == desc:
                x[key] = x['ROMol'].apply(func)
    return x


def normalize_df(df, bound, scale=5) : 
    for prop, b in bound.items() : 
        df = df.astype({prop : 'float64'})
        min_val, max_val = df[prop].min(), df[prop].max()
        df[prop] = (df[prop] - min_val) / (max_val - min_val) if max_val > min_val else 0
        bound[prop][0] = ((b[0] - min_val) / (max_val - min_val)) * scale + scale if max_val > min_val else 0
        bound[prop][1] = ((b[1] - min_val) / (max_val - min_val)) * scale + scale if max_val > min_val else 0
    bound = np.array(list(bound.values()))
    lower_bound, upper_bound = list(bound[:,0]), list(bound[:,1])
    lower_bound.append(lower_bound[0])
    upper_bound.append(upper_bound[0])
    return df, lower_bound, upper_bound


def get_df_mean_std(df, scale=5) : 
    stats = df.describe().T
    stats['mean'] = stats['mean'] * scale + scale
    stats['std'] = stats['std'] * scale
    return stats[['mean','std']]




def define_radial_axes_angles(n_axes):
    """Define angles (radians) for radial (x-)axes."""
    angles = [i / float(n_axes) * 2 * math.pi for i in range(n_axes)]
    return angles + angles[:1]

def bound_to_label(bound):
    """Convert bounds dictionary to formatted labels."""
    return [f'{k} {v[0], v[1]}' for k, v in bound.items()]




def molecular_radar_plot(dataframe, bound, scale=5, output_path=None) : 
    properties_labels = bound_to_label(bound)
    angles = define_radial_axes_angles(len(bound))

    norm_df, lower_bound, upper_bound = normalize_df(dataframe, bound, scale)
    ymax = int(max(upper_bound)) + 1

    stats_df = get_df_mean_std(norm_df, scale)
    stats_df = pd.concat([stats_df, stats_df.head(1)])

    plt.figure(figsize=(6,6)) 
    ax = plt.subplot(111, polar=True)
    ax.fill_between(angles, lower_bound, upper_bound, color='cornflowerblue', alpha=0.2)
    ax.plot(angles, stats_df["mean"], "b", lw=3, ls="-")
    ax.plot(angles, stats_df["mean"] + stats_df["std"], "orange", lw=1, ls="--")
    ax.plot(angles, stats_df["mean"] - stats_df["std"], "orange", lw=1, ls="--")
    ax.legend(('boundary', "mean", "mean ± std"), loc=(1.3, 0.7), labelspacing=0.3, fontsize=16)
    ax.set_theta_offset(math.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_rlabel_position(180)
    plt.xticks(angles, [])
    plt.ylim(0, ymax)
    plt.yticks(range(1, ymax), [])


    for _, (angle, label) in enumerate(zip(angles[:-1], properties_labels)):
        if angle == 0: ha = "center"
        elif 0 < angle < math.pi: ha = "left"
        elif angle == math.pi: ha = "center"
        else: ha = "right"
        ax.text(
            x=angle,
            y=ymax + 1,
            s=label,
            size=16, 
            horizontalalignment=ha,
            verticalalignment="center")
        
    if output_path :
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        







# Parsing arguments function
def parse_exist(x) : 
    dic = {}
    if x is not None :
        x = x.split(' ')
        for i in x : 
            desc, bound = i.split(':')
            lower, upper = bound[1:-1].split(',')
            dic[desc] = [float(lower), float(upper)]
    return dic

def parse_bound(x) : 
    min, max = x.split(',')
    return [float(min), float(max)]
