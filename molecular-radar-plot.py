import argparse
import rdkit 
from rdkit import Chem 
from rdkit.Chem import Descriptors, PandasTools
from utils import * 
# MW:[200,300] HBA:[2,5] HBD:[1,3]

parser = argparse.ArgumentParser(description='Molecular Radar Plot')
parser.add_argument('-i', '--input', type=str, help='Path to the file containing SMILES')
parser.add_argument('-e', '--exist', type=parse_exist, help='Descriptors that already exist in the .sdf file. Example: "MW:[200,300] HBA:[2,5] HBD:[1,3]"')
parser.add_argument('-o', '--output', type=str, help='Path to the output file (.png)')
for d, _ in Descriptors._descList : 
    parser.add_argument(f'--{d}', type=parse_bound)
args = parser.parse_args()



input = args.input
bound = vars(args) 
del bound['input']
del bound['exist']
for key in list(bound.keys()):
    if bound[key] is None:
        del bound[key]


try : 
    dataframe = PandasTools.LoadSDF(input, molColName='ROMol')
    if args.exist : 
        for key in args.exist.keys() : 
            if key not in dataframe.columns :
                print(f'{key} not in the dataframe'); exit()
    if bound : 
        dataframe = update_dataframe(dataframe, bound)
except: 
    smiles_list = read_smi(input)
    dataframe = update_dataframe(smiles_list, bound)

bound.update(args.exist)



molecular_radar_plot(dataframe, bound, output_path=args.output if args.output else 'molecular_radar_plot.png')