import os
import csv
import glob
from tqdm import tqdm
import pandas as pd

# Set the directory containing your CSV files
directory = '/home/user/LabData/Reza/data/Wout/121024/res/'
# directory = '/home/user/LabData/Reza/data/Wout/120924/nf_output/res/'

# Name of the output file
output_file = '/home/user/LabData/Reza/data/Wout/121024/combined.csv'
# output_file = '/home/user/LabData/Reza/data/Wout/121024/res/combined_old.csv'

# Get a list of all CSV files in the directory
all_files = glob.glob(os.path.join(directory, "*.csv"))

# Process the first file separately to get the header
first_file = all_files[0]
df = pd.read_csv(first_file)

columns = df.columns
df = pd.DataFrame(columns=columns)
df.to_csv(output_file, index=False)
# header = True
# Process the rest of the files 
for file in tqdm(all_files):
    df = pd.read_csv(file)
    df.to_csv(output_file, mode='a', header=False, index=False)
    # header = False

print(f"Merged {len(all_files)} CSV files into {output_file}")