import os
import re
import csv
from tqdm import tqdm
import argparse


def get_file_size(file_path):
    """Get the file size in bytes."""
    try:
        return os.path.getsize(file_path)
    except:
        return 0


def extract_values(file_path):
    """Extract decoy and ddg_norepack values from the text file with progress display."""
    decoy_values = []
    ddg_norepack_values = []
    error_lines = []
    decoy_pattern = re.compile(r'"decoy"\s*:\s*"([^"]+)"')
    ddg_norepack_pattern = re.compile(r'"ddg_norepack"\s*:\s*([-+]?\d*\.\d+|[-+]?\d+)')

    file_size = get_file_size(file_path)
    processed_size = 0

    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            # Get the total number of lines for progress display
            total_lines = sum(1 for _ in file)
            file.seek(0)  # Reset the file pointer

            # Create a progress bar
            with tqdm(total=total_lines, desc="Processing file", unit="line") as pbar:
                for line_num, line in enumerate(file, 1):
                    line = line.strip()
                    if not line:
                        pbar.update(1)
                        continue

                    # Use regular expressions to find decoy and ddg_norepack values
                    decoy_matches = decoy_pattern.findall(line)
                    ddg_norepack_matches = ddg_norepack_pattern.findall(line)

                    if decoy_matches and ddg_norepack_matches:
                        decoy_values.extend(decoy_matches)
                        ddg_norepack_values.extend(ddg_norepack_matches)

                    # Update the progress bar
                    processed_size += len(line.encode('utf-8')) + 1  # Add the newline character
                    pbar.set_postfix({"Found Decoys": len(decoy_values)})
                    pbar.update(1)
    except Exception as e:
        print(f"Error processing file: {e}")
        return [], [], []

    return decoy_values, ddg_norepack_values, error_lines


def write_to_csv(decoy_values, ddg_norepack_values, output_file, append=False):
    """Write decoy and ddg_norepack values to a CSV file."""
    mode = 'a' if append else 'w'
    with open(output_file, mode, encoding='utf-8', newline='') as csvfile:
        fieldnames = ['decoy', 'ddg_norepack']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write the header only if the file is opened in write mode
        if mode == 'w':
            writer.writeheader()

        # Write the data rows
        for decoy, ddg_norepack in zip(decoy_values, ddg_norepack_values):
            writer.writerow({'decoy': decoy, 'ddg_norepack': ddg_norepack})


def main():
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Extract decoy and ddg_norepack values from a text file.')
    parser.add_argument('-i', '--input', type=str, required=True, help='Path to the input file.')
    parser.add_argument('-o', '--output', type=str, required=True, help='Path to the output CSV file.')
    parser.add_argument('-e', '--exist_file', action='store_true',
                        help='If specified, assume the output file already exists and append data to it.')

    # Parse the arguments
    args = parser.parse_args()
    file_path = args.input
    output_file = args.output
    append = args.exist_file

    # Check if the input file exists
    if not os.path.exists(file_path):
        print(f"Error: The file {file_path} does not exist.")
        return

    print(f"Starting to extract decoy and ddg_norepack values from {file_path}...")

    # Extract values (with progress bar)
    decoy_values, ddg_norepack_values, error_lines = extract_values(file_path)

    # Write to the output CSV file
    write_to_csv(decoy_values, ddg_norepack_values, output_file, append)

    print(f"Successfully extracted {len(decoy_values)} decoy and ddg_norepack pairs and wrote them to {output_file}")
    if error_lines:
        print(f"Note: There are {len(error_lines)} lines containing 'decoy' or 'ddg_norepack' but with a mismatched format.")


if __name__ == "__main__":
    main()