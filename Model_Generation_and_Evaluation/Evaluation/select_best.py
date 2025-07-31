import csv
import argparse

# Store the row with the minimum ddg_norepack for each identifier
best_rows = {}

# Parse command line arguments
parser = argparse.ArgumentParser(description='Select rows with the minimum ddg_norepack based on specific identifiers.')
parser.add_argument('-i', '--input', required=True, help='Input CSV file path')
parser.add_argument('-o', '--output', required=True, help='Output CSV file path')
args = parser.parse_args()

# Read the input CSV file
with open(args.input, 'r', newline='') as csvfile:
    reader = csv.reader(csvfile)
    header = next(reader)  # Read header row
    for row in reader:
        decoy = row[0]
        # Extract characters 6-9 (indexes 5-8) as the identifier
        identifier = decoy[5:9] if len(decoy) >= 9 else decoy
        # Check if ddg_norepack value exists
        if row[1] != '':
            ddg_norepack = float(row[1])
        else:
            continue  # Skip rows with empty ddg_norepack
        
        # Update best row for the current identifier
        if identifier not in best_rows:
            best_rows[identifier] = row
        else:
            # Get current best ddg value
            current_best_ddg = float(best_rows[identifier][1]) if best_rows[identifier][1] != '' else float('inf')
            if ddg_norepack < current_best_ddg:
                best_rows[identifier] = row

# Write the selected rows to the output CSV file
with open(args.output, 'w', newline='') as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(header)  # Write header
    for row in best_rows.values():
        writer.writerow(row)

print(f"Filtering completed. Results saved to {args.output}")
