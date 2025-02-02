# Molecular Radar Plot

A flexible Python tool for generating radar plots to visualize molecular properties using RDKit descriptors. This tool allows you to visualize and analyze the distribution of any molecular descriptor available in RDKit across your chemical dataset.
---

## Key Features

- Support for multiple input formats (SMILES, SDF)
- Visualization of molecular property distributions with mean and standard deviation
- Flexible descriptor selection using any RDKit molecular descriptor
- Customizable property ranges for visualization
- Support for pre-calculated descriptors in SDF files
---

## Installation

### Prerequisites
```bash
pip install rdkit pandas matplotlib numpy tqdm
```
---

## Usage

### Command Line Arguments

- `-i, --input`: Path to input file (supports .smi, .txt, or .sdf)
- `-o, --output`: Path to output plot (default: molecular_radar_plot.png)
- `-e, --exist`: Descriptors that already exist in the SDF file (format: "DESC1:[min1,max1] DESC2:[min2,max2]")
- `--{descriptor}`: Any RDKit descriptor with its range (format: min,max)

To view available rdkit Descriptors:
```bash
from rdkit.Chem import Descriptors
print(Descriptors._descList)
```

Some common examples:
- `ExactMolWt`: Molecular Weight
- `TPSA`: Topological Polar Surface Area
- `MolLogP`: Calculated LogP
- `NumHAcceptors`: Number of H-Bond Acceptors
- `NumHDonors`: Number of H-Bond Donors
- `NumRotatableBonds`: Number of Rotatable Bonds
- `qed`: Quantitative estimate of drug-likeness

Basic usage with descriptors available for calculation in rdkit:
```bash
python molecular-radar-plot.py -i molecules.smi --ExactMolWt 200,300 --NumHAcceptors 2,5 --NumHDonors 1,3 --qed 0.7,0.9
```

Using SDF file with existing descriptors:
```bash
python molecular-radar-plot.py -i molecules.sdf -e "MW:[200,500] HBA:[2,5] HBD:[1,3]"
```

---

## Output

The tool generates a radar plot showing:
- Blue shaded area: Property boundaries
- Blue solid line: Mean values
- Orange dashed lines: Standard deviation bounds

![Example Radar Plot](data/radar-plot.png)

---


## Contributing

Feel free to open issues or submit pull requests with improvements or bug fixes.
