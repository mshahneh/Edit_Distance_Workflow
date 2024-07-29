import os
import csv
from tqdm import tqdm

# Set the directory containing your CSV files
directory = '/home/user/LabData/Reza/data/Wout/nf_output/res'

# Name of the output file
output_file = '/home/user/LabData/Reza/data/Wout/combined_small.csv'

# Get a list of CSV files in the directory
csv_files = [f for f in os.listdir(directory) if f.endswith('.csv')]
print(len(csv_files), "files found.")

# Initialize a variable to track whether the header has been written
header_written = False

with open(output_file, 'w', newline='') as outfile:
    writer = csv.writer(outfile)
    
    for filename in tqdm(csv_files):
        with open(os.path.join(directory, filename), 'r') as infile:
            reader = csv.reader(infile)
            header = next(reader)  # Read the header
            
            if not header_written:
                writer.writerow(header)  # Write the header only once
                header_written = True
                
            for row in reader:
                writer.writerow(row)  # Write the data

print("Files have been combined into:", output_file)