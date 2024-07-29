
#from ipywidgets import interact, fixed, widgets
import pandas as pd
from rdkit import Chem, Geometry
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
import argparse
from tqdm import tqdm
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
import os
from modifinder.utilities import mol_utils as mu

def solve_pair(smiles1, smiles2):
    if type(smiles1) == str:
        mol1 = Chem.MolFromSmiles(smiles1)
    else:
        mol1 = smiles1
    
    if type(smiles2) == str:
        mol2 = Chem.MolFromSmiles(smiles2)
    else:
        mol2 = smiles2
    
    if mol1.GetNumAtoms() > 60 or mol2.GetNumAtoms() > 60:
        raise ValueError("The molecules are too large.")
    
    fpgen = AllChem.GetRDKitFPGenerator(maxPath=3,fpSize=512)
    fps = [fpgen.GetFingerprint(x) for x in [mol1, mol2]]
    tanimoto = DataStructs.TanimotoSimilarity(fps[0],fps[1])
    if tanimoto < 0.2:
        raise ValueError("The Tanimoto similarity is too low.")
    
    mcs1 = rdFMCS.FindMCS([mol1, mol2])
    mcs_mol = Chem.MolFromSmarts(mcs1.smartsString)
    if mcs_mol.GetNumAtoms() < mol1.GetNumAtoms()//2 and mcs_mol.GetNumAtoms() < mol2.GetNumAtoms()//2:
        raise ValueError("The MCS is too small.")
    if mcs_mol.GetNumAtoms() <= 2:
        raise ValueError("The MCS is too small.")
    
    dist1 = mu.get_modification_edges(mol1, mcs_mol)
    dist2 = mu.get_modification_edges(mol2, mcs_mol)
    distance = len(dist1) + len(dist2)

    is_sub = ((len(dist1) == 0) or (len(dist2) == 0))
    return distance, tanimoto, is_sub

def main(data_x, out_dir, data_y):
    mols = {i: Chem.MolFromSmiles(row['Smiles']) for i, row in data_y.iterrows()}

    count = 0
    with open(out_dir, 'w') as f:
        f.write('smiles1,smiles2,distance,tanimoto,is_sub,inchi1,inchi2\n')
        for i, row1 in tqdm(data_x.iterrows(), total=len(data_x), ascii=True):
            moli = Chem.MolFromSmiles(row1['Smiles'])
            if moli.GetNumAtoms() > 60:
                continue
            for j, row2 in data_y.iterrows():
                # only for upper triangular matrix
                if i > j:
                    continue
                try:
                    distance, tanimoto, is_sub = solve_pair(moli, mols[j])
                    f.write(f"{row1['Smiles']},{row2['Smiles']},{distance},{tanimoto},{is_sub},\"{row1['INCHI']}\",\"{row2['INCHI']}\"\n")
                    # print("here")
                except Exception as e:
                    # print(e)
                    continue
        
        count += 1
        # if 10 percent of progress is made, save the file
        if count % (len(data_x) // 10) == 0:
            f.flush()
            os.fsync(f.fileno())

def get_data_index(data, index, batch_count):
    batch_size = len(data) // batch_count
    if index < len(data) % batch_count:
        batch_size += 1
        start = index * batch_size
        end = start + batch_size
    else:
        start = index * batch_size + len(data) % batch_count
        end = start + batch_size
    
    if end > len(data):
        end = len(data)
    
    return start, end

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate edit distance and Tanimoto similarity between pairs of molecules.')
    parser.add_argument('data', type=str, help='Path to the input data file.')
    parser.add_argument('out_dir', type=str, help='Path to the output directory.')
    parser.add_argument('--batch_x', type=int, help='Number of rows in the batch.', default=1)
    parser.add_argument('--batch_y', type=int, help='Number of columns in the batch.', default=1)
    parser.add_argument('--batch_index', type=int, help='Index of the batch to process.', default=0)
    parser.add_argument('--cached_dir', type=str, help='Path to the directory containing cached data.', default=None)
    args = parser.parse_args()

    # if output dir is relative, make it absolute
    if not os.path.isabs(args.out_dir):
        args.out_dir = os.path.join(os.getcwd(), args.out_dir)
    
    if not os.path.exists(os.path.dirname(args.out_dir)):
        os.makedirs(os.path.dirname(args.out_dir))
    
    if args.cached_dir is not None:
        print("Checking cache...")
        file_name = os.path.basename(args.out_dir)
        print("file_name", file_name)
        cached_file = os.path.join(args.cached_dir, file_name)
        print(cached_file)
        if os.path.exists(cached_file):
            print("Cache found!.")
            # rewrite the file in the output directory
            print(args.out_dir, "is being written.")
            with open(cached_file, 'r') as f:
                with open(args.out_dir, 'w') as f2:
                    f2.write(f.read())
            exit(0)
    else:
        print("No cache directory provided.")
    

    data = pd.read_csv(args.data)
    batch_x = args.batch_x
    batch_y = args.batch_y
    batch_index = args.batch_index
    out_dir = args.out_dir

    index_x = batch_index // batch_y
    index_y = batch_index % batch_y

    data_y_start, data_y_end = get_data_index(data, index_y, batch_y)
    data_start, data_end = get_data_index(data, index_x, batch_x)

    # only do for upper triangular matrix
    data_y_start = max(data_y_start, data_start)

    data_y = data.iloc[data_y_start:data_y_end]
    data_x = data.iloc[data_start:data_end]
    del data

    main(data_x, out_dir, data_y)


