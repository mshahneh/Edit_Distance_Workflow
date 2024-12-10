""" 
the cached mols are a dictionary with inchi as key and rdkit mol object as value
the cached motifs is a list of motifs molecules
""" 

import pandas as pd
from rdkit import Chem
import pickle
import argparse
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create pickled cache files for the input data.')
    parser.add_argument('mols', type=str, help='Path to the input data file.')
    parser.add_argument('motifs', type=str, help='Path to the motifs file.')
    parser.add_argument('out_dir', type=str, help='Path to the output directory.')
    
    args = parser.parse_args()
    
    if not os.path.isabs(args.out_dir):
        args.out_dir = os.path.join(os.getcwd(), args.out_dir)
        
    if not os.path.exists(os.path.dirname(args.out_dir)):
        os.makedirs(os.path.dirname(args.out_dir))
        
    mols = pd.read_csv(args.mols)
    motifs = pd.read_csv(args.motifs)
    
    cached_mols = dict()
    for index, row in mols.iterrows():
        try:
            cached_mols[row['INCHI']] = Chem.MolFromSmiles(row['Smiles'])
        except:
            cached_mols[row['INCHI']] = None
    
    with open(args.out_dir + '/mols.pkl', 'wb') as f:
        pickle.dump(cached_mols, f)
    
    smart_mols = motifs['smarts'].apply(Chem.MolFromSmarts)
    smart_mols = list(smart_mols)
    with open(args.out_dir + '/motifs_mols.pkl', 'wb') as f:
        pickle.dump(smart_mols, f)
    
        