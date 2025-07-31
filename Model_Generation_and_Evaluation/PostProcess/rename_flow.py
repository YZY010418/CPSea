import os
import shutil
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Process and rename PDB files from a specified folder.')
parser.add_argument('-i', '--input', required=True, help='Path to the main folder containing subfolders with PDB files')
args = parser.parse_args()

# Define main folder path and target folder paths
main_folder = args.input
flow_gt_folder = "flow_gt"
flow_renamed_folder = "flow_renamed"

# Create target folders if they don't exist
os.makedirs(flow_gt_folder, exist_ok=True)
os.makedirs(flow_renamed_folder, exist_ok=True)

# Iterate through all subfolders in the main folder
for sub_folder in os.listdir(main_folder):
    sub_folder_path = os.path.join(main_folder, sub_folder)
    # Check if it's a directory
    if os.path.isdir(sub_folder_path):
        # Iterate through all files in the subfolder
        for file_name in os.listdir(sub_folder_path):
            if file_name.endswith(".pdb"):
                old_file_path = os.path.join(sub_folder_path, file_name)
                if file_name == "gt.pdb":
                    # Process gt.pdb files
                    new_file_name = f"{sub_folder}_gt.pdb"
                    new_file_path = os.path.join(flow_gt_folder, new_file_name)
                else:
                    # Process numbered PDB files (extract numeric part)
                    # Assumes filename format contains numbers (e.g., xxx123.pdb, 123.pdb, etc.)
                    # Extract numbers from filename as the identifier
                    num = ''.join(filter(str.isdigit, file_name))
                    if not num:
                        print(f"File {file_name} has no number identifier, skipping processing")
                        continue
                    new_file_name = f"flow_{sub_folder}_{num}.pdb"
                    new_file_path = os.path.join(flow_renamed_folder, new_file_name)
                
                # Copy and rename the file instead of moving
                try:
                    shutil.copy2(old_file_path, new_file_path)  # Using copy2 to preserve metadata
                    print(f"Copied and renamed {old_file_path} to {new_file_path}")
                except Exception as e:
                    print(f"Error processing {old_file_path}: {e}")
    
