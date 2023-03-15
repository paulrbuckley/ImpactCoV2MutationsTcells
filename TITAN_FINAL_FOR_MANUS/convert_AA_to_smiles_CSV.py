import sys
import csv
import os
import itertools
from pytoda.proteins import aas_to_smiles

# Extract command line arguments, remove file extension and attach to output_filename
input_filename1 = sys.argv[1]
input_filename2 = os.path.splitext(input_filename1)[0]
filenames = (input_filename2, "SMILES_OUT.csv")
output_filename = "".join(filenames)

def main(argv):
    with open(input_filename1, 'r') as f:
        with open(output_filename, 'w') as g:
    
            csvread = csv.reader(f)
            print(csvread)
            csvwrite = csv.writer(g)
    
            header = next(csvread)
            header.append("SMILES")
            csvwrite.writerow(header)
  
            for row in csvread:
                peptide = row[0]
                smiles = aas_to_smiles(peptide)
                row.append(smiles)
                csvwrite.writerow(row)

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
