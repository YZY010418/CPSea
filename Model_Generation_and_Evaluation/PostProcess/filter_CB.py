import os
import shutil
import csv
import sys

def move_files_based_on_csv(csv_path, source_dir, target_dir):
    # Ensure the target directory exists
    os.makedirs(target_dir, exist_ok=True)
    
    # Read the CSV file and get the filenames to be moved
    files_to_move = set()
    with open(csv_path, 'r', encoding='utf-8') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if len(row) >= 2:
                filename = row[0]
                try:
                    number = float(row[1])
                    if not (3 <= number <= 8):
                        files_to_move.add(filename)
                except ValueError:
                    print(f"Warning: Unable to convert value '{row[1]}' to a number, skipping this row.")
    
    # Traverse the source directory and its subdirectories to move matching files
    moved_count = 0
    for root, _, files in os.walk(source_dir):
        for file in files:
            if file in files_to_move:
                source_path = os.path.join(root, file)
                target_path = os.path.join(target_dir, file)
                
                # Handle the case where the target file already exists
                counter = 1
                original_name, ext = os.path.splitext(target_path)
                while os.path.exists(target_path):
                    target_path = f"{original_name}_{counter}{ext}"
                    counter += 1
                
                try:
                    shutil.move(source_path, target_path)
                    moved_count += 1
                    print(f"Moved: {source_path} -> {target_path}")
                except Exception as e:
                    print(f"Error: Failed to move file {source_path}: {e}")
    
    print(f"Operation completed. Total {moved_count} files moved.")

if __name__ == "__main__":
    # Check if enough arguments are provided
    if len(sys.argv) != 4:
        print("Usage: python filter_CB.py <csv_file_path> <source_directory> <target_directory>")
        sys.exit(1)
    
    # Get parameters from command line arguments
    csv_path = sys.argv[1]
    source_dir = sys.argv[2]
    target_dir = sys.argv[3]
    
    move_files_based_on_csv(csv_path, source_dir, target_dir)
