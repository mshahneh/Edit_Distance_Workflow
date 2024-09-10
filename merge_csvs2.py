import os
import csv
import glob
from tqdm import tqdm
import pandas as pd

# Set the directory containing your CSV files
directory = '/home/user/LabData/Reza/data/Wout/res_fixed'

# Name of the output file
output_file = '/home/user/LabData/Reza/data/Wout/combined2.csv'

# Get a list of all CSV files in the directory
all_files = glob.glob(os.path.join(directory, "*.csv"))

# Process the first file separately to get the header
first_file = all_files[0]


# Write the header to the output file
columns = ["smiles1","smiles2","distance","tanimoto","is_sub","inchi1","inchi2"]
df = pd.DataFrame(columns=columns)
df.to_csv(output_file, index=False)

# Process the rest of the files 
for file in tqdm(all_files):
    # read csv file that doesnt have header
    df = pd.read_csv(file, header=None)
    df.to_csv(output_file, mode='a', header=False, index=False)

print(f"Merged {len(all_files)} CSV files into {output_file}")