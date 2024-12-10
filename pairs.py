
#from ipywidgets import interact, fixed, widgets
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
import argparse
from tqdm import tqdm
from rdkit import DataStructs
import os
import mol_utils as mu
from myopic_mces import MCES
from motif_base_edit_distance import get_min_motif_to_cover
import pickle

atom_num_threshold = 55
upper_limit_edit_distance = 10


def motif_edit_distance(mol, mcs_mol, motifs, top_limit=None):
    # print('here and the limit is', top_limit)
    modifications = mu.get_modification_graph(mol, mcs_mol)
    # print("modifications", modifications)
    distance = 0
    for modification in modifications:
        res, sections = get_min_motif_to_cover(modification[0], motifs, top_limit)
        distance += max(res, len(modification[1]))
        if distance > upper_limit_edit_distance:
            return distance
    # print("done with motif edit distance")
    return distance
        

def solve_pair(smiles1, smiles2, motifs_dict, top_limit=None, important_columns=None):
    if isinstance(smiles1, str):
        mol1 = Chem.MolFromSmiles(smiles1)
    else:
        mol1 = smiles1
    
    if isinstance(smiles2, str):
        mol2 = Chem.MolFromSmiles(smiles2)
    else:
        mol2 = smiles2
    
    if mol1.GetNumAtoms() > atom_num_threshold or mol2.GetNumAtoms() > atom_num_threshold:
        raise ValueError("The molecules are too large.")
    
    # print("going to calculate tanimoto")
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024) for mol in [mol1, mol2]]
    tanimoto = DataStructs.TanimotoSimilarity(fps[0],fps[1])
    if tanimoto < 0.1:
        raise ValueError("The Tanimoto similarity is too low.")
    
    # print("tanimoto", tanimoto)
    # print("going to calculate mcs")
    
    mcs1 = rdFMCS.FindMCS([mol1, mol2])
    mcs_mol = Chem.MolFromSmarts(mcs1.smartsString)
    
    # print("mcs", mcs_mol.GetNumAtoms())
    if mcs_mol.GetNumAtoms() < mol1.GetNumAtoms()//2 or mcs_mol.GetNumAtoms() < mol2.GetNumAtoms()//2:
        raise ValueError("The MCS is too small.")
    if mcs_mol.GetNumAtoms() <= 2:
        raise ValueError("The MCS is too small.")
    
    # print("going to calculate edit distance")
    dist1, dist2 = mu.get_edit_distance_detailed(mol1, mol2, mcs_mol)
    
    distance = dist1 + dist2
    is_sub = ((dist1 == 0) or (dist2 == 0))
    # print("distance", distance)
    

    res = dict()
    res['distance'] = distance
    res['tanimoto'] = tanimoto
    res['is_sub'] = is_sub
    res['added'] = dist1
    res['removed'] = dist2
    
    if 'mces' in important_columns:
        if not isinstance(smiles1, str):
            smiles1 = Chem.MolToSmiles(mol1)
        if not isinstance(smiles2, str):
            smiles2 = Chem.MolToSmiles(mol2)
        mces_value = MCES(smiles1, smiles2, threshold=upper_limit_edit_distance)
        res['mces'] = round(mces_value[1], 2)
    
    if "motif_distance" in important_columns:
        # print("going to calculate motif edit distance")
        if distance >= upper_limit_edit_distance:
            motif_removed = dist2
            motif_added = dist1
        else:
            # print("not passed the threshold")
            motif_removed = motif_edit_distance(mol1, mcs_mol, motifs_dict, top_limit)
            motif_added = motif_edit_distance(mol2, mcs_mol, motifs_dict, top_limit)
        res['motif_distance'] = min(upper_limit_edit_distance, max(distance, motif_removed + motif_added))
        res['motif_removed'] = motif_removed
        res['motif_added'] = motif_added
    
    res['distance'] = min(upper_limit_edit_distance, distance)
    
    return res

def main(cached_structures, cached_motifs, top_limit, data_x, out_dir, data_y, verbose=True):
    important_columns = ['distance','tanimoto','is_sub','added','removed','mces','motif_distance','motif_removed','motif_added']

    with open(out_dir, 'w') as f:
        # f.write('distance,tanimoto,is_sub,inchi1,inchi2\n')
        f.write("inchi1,inchi2," + ','.join(important_columns) + '\n')
        for i, row1 in tqdm(data_x.iterrows(), total=len(data_x), ascii=True, disable=not verbose):
            ###### DEBUGGING ######
            # if i != 76:
            #     continue
            
            # print(i)
            moli = cached_structures[row1['INCHI']]
            if moli.GetNumAtoms() > atom_num_threshold:
                continue
            for j, row2 in data_y.iterrows():
                ###### DEBUGGING ######
                # if j != 1348:
                #     continue
                
                # print(row1['Smiles'], row2['Smiles'])
                
                molj = cached_structures[row2['INCHI']]
                if molj.GetNumAtoms() > atom_num_threshold:
                    continue
                # warning!!! this is to make sure the process is only called for the upper triangular matrix, this assumes
                # that the data is not reindexed and the data_x and data_y are just slices of the original dataframe
                if i > j:
                    continue
                try:
                    res = solve_pair(moli, molj, cached_motifs, top_limit, important_columns)
                    str_res = [str(res[col]) for col in important_columns]
                    f.write(f"\"{row1['INCHI']}\",\"{row2['INCHI']}\",")
                    f.write(','.join(str_res) + '\n')
                    # f.write(f"{distance},{tanimoto},{is_sub},\"{row1['INCHI']}\",\"{row2['INCHI']}\"\n")
                    # f.write(f"{distance},{tanimoto},{is_sub},\"{row1['INCHI']}\",\"{row2['INCHI']}\",{mces},{res['added']},{res['removed']}\n")
                    # print("here")
                except Exception as e:
                    # f.write("Error: " + str(e) + '\n')
                    if 'int' in e.args[0] or 'float' in e.args[0]:
                        raise e
                    continue
            
            f.flush()
        
        
        # count += 1
        # # if 10 percent of progress is made, save the file
        # if count in save_intervals:
        #     f.flush()
        #     os.fsync(f.fileno())

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
    parser.add_argument('cached_structures', type=str, help='Path to the cached structures file.')
    parser.add_argument('motifs', type=str, help='Path to the motifs file.')
    parser.add_argument('out_dir', type=str, help='Path to the output directory.')
    parser.add_argument('--batch_x', type=int, help='Number of rows in the batch.', default=1)
    parser.add_argument('--batch_y', type=int, help='Number of columns in the batch.', default=1)
    parser.add_argument('--batch_index', type=int, help='Index of the batch to process.', default=0)
    parser.add_argument('--cached_dir', type=str, help='Path to the directory containing cached data.', default=None)
    parser.add_argument('--top_limit', type=int, help='limit for depth of motif edit distance', default=4)
    parser.add_argument('--verbose', action='store_true', help='Print verbose output.')
    args = parser.parse_args()

    # if output dir is relative, make it absolute
    if not os.path.isabs(args.out_dir):
        args.out_dir = os.path.join(os.getcwd(), args.out_dir)
        
    print(args.out_dir)
    
    if not os.path.exists(os.path.dirname(args.out_dir)):
        os.makedirs(os.path.dirname(args.out_dir))
    
    if args.cached_dir is not None:
        # print("Checking cache...")
        file_name = os.path.basename(args.out_dir)
        # print("file_name", file_name)
        cached_file = os.path.join(args.cached_dir, file_name)
        # print(cached_file)
        if os.path.exists(cached_file):
            if args.verbose:
                print("Cache found!.")
            # rewrite the file in the output directory
            # print(args.out_dir, "is being written.")
            with open(cached_file, 'r') as f:
                with open(args.out_dir, 'w') as f2:
                    f2.write(f.read())
            exit(0)
    else:
        if args.verbose:    
            print("No cache directory provided.")
    
    with open(args.cached_structures, 'rb') as f:
        cached_structures = pickle.load(f)
        
    with open(args.motifs, 'rb') as f:
        motifs_dict = pickle.load(f)
    

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

    if args.verbose:
        print("done reading data", data_x.shape, data_y.shape)
    main(cached_structures, motifs_dict, args.top_limit, data_x, out_dir, data_y, args.verbose)

