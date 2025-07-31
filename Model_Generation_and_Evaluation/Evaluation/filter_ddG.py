import argparse
import csv
import os
import shutil

def main():
    parser = argparse.ArgumentParser(description='Move corresponding PDB files based on values in CSV file')
    parser.add_argument('csv_file', help='Path to the CSV file')
    parser.add_argument('input_folder', help='Path to the input folder containing PDB files')
    parser.add_argument('output_folder', help='Path to the output folder for moved files')
    
    args = parser.parse_args()
    
    os.makedirs(args.output_folder, exist_ok=True)
    
    pdb_basenames = set()
    for filename in os.listdir(args.input_folder):
        if filename.endswith('.pdb'):
            basename = os.path.splitext(filename)[0]
            pdb_basenames.add(basename)
    
    with open(args.csv_file, 'r') as csvfile:
        csvreader = csv.reader(csvfile)
        
        for row_num, row in enumerate(csvreader, 1):
            if len(row) < 2:
                print(f"Warning: Insufficient data in row {row_num}, skipped")
                continue
            
            first_column = row[0]
            if '_0001' in first_column:
                identifier = first_column.split('_0001')[0]
            else:
                print(f"Warning: First column in row {row_num} does not contain '_0001', skipped")
                continue
            
            if identifier not in pdb_basenames:
                print(f"Warning: No corresponding PDB file for identifier '{identifier}', row {row_num} skipped")
                continue
            
            try:
                value = float(row[1])
            except ValueError:
                print(f"Warning: Second column in row {row_num} is not a valid float, skipped")
                continue
            
            if value >= 0:
                pdb_filename = f"{identifier}.pdb"
                source_path = os.path.join(args.input_folder, pdb_filename)
                dest_path = os.path.join(args.output_folder, pdb_filename)
                
                shutil.move(source_path, dest_path)
                print(f"Moved file: {pdb_filename}")
    
    print("Processing completed!")

if __name__ == "__main__":
    main()