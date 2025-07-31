import sys
import csv
import math
import numpy as np
from pathlib import Path
from scipy import stats

def analyze_interaction_data(file_path):
    # Check if the file exists
    file = Path(file_path)
    if not file.exists():
        print(f"Error: The file '{file_path}' does not exist!")
        return 1
    
    # Check if the file extension is CSV
    if file.suffix.lower() != '.csv':
        print(f"Error: Please provide a CSV file, not '{file.suffix}'")
        return 1
    
    # Define interaction type columns
    interaction_types = [
        'hydrophobic_interactions', 'hydrogen_bonds', 'water_bridges',
        'salt_bridges', 'pi_stacks', 'pi_cation_interactions',
        'halogen_bonds', 'metal_complexes'
    ]
    
    try:
        # Initialize data structures
        interaction_data = {col: [] for col in interaction_types}
        ratio_data = {col: [] for col in interaction_types}
        
        # Read the CSV file and process the data
        print(f"Reading file: {file_path}")
        print("Processing data...")
        with open(file_path, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            
            # Check if all required columns exist
            if not all(col in reader.fieldnames for col in interaction_types):
                missing = [col for col in interaction_types if col not in reader.fieldnames]
                print(f"Error: Missing necessary columns: {', '.join(missing)}")
                return 1
            
            # Process each row
            for i, row in enumerate(reader):
                # Extract interaction data
                values = [int(row.get(col, 0)) for col in interaction_types]
                total = sum(values)
                
                # Save raw data
                for col_idx, col in enumerate(interaction_types):
                    interaction_data[col].append(values[col_idx])
                
                # Calculate and save ratios
                if total > 0:
                    for col_idx, col in enumerate(interaction_types):
                        ratio = values[col_idx] / total
                        ratio_data[col].append(ratio)
                
                # Display progress every 1 million rows
                if (i + 1) % 1000000 == 0:
                    print(f"{i + 1} rows processed")
        
        # Calculate statistical data
        print("Calculating statistical data...")
        stats_data = []
        for col in interaction_types:
            ratios = ratio_data[col]
            if not ratios:
                mean_val = 0
                std_val = 0
            else:
                mean_val = sum(ratios) / len(ratios)
                variance = sum((x - mean_val) ** 2 for x in ratios) / (len(ratios) - 1) if len(ratios) > 1 else 0
                std_val = math.sqrt(variance)
            
            stats_data.append({
                'Interaction_Type': col,
                'Mean': mean_val,
                'StdDev': std_val
            })
        
        # Print statistical results
        print("\nStatistical Results:")
        print(f"{'Interaction Type':<25} {'Mean':<10} {'Std Dev':<10}")
        print("-" * 50)
        for data in stats_data:
            print(f"{data['Interaction_Type']:<25} {data['Mean']:10.4f} {data['StdDev']:10.4f}")
        
        # Save statistical results to CSV
        stats_file = file.with_name(f"{file.stem}_stats.csv")
        write_dict_to_csv(stats_file, stats_data)
        print(f"\nStatistical results saved to: {stats_file}")
        
        # Calculate and save KDE data
        print("\nCalculating KDE data...")
        kde_data = generate_kde_data(ratio_data)
        kde_file = Path("CPSet_plip_KDE.csv")
        write_dict_to_csv(kde_file, kde_data)
        print(f"KDE data saved to: {kde_file}")
        
        return 0
    
    except Exception as e:
        print(f"An error occurred while processing the file: {str(e)}")
        return 1

def write_dict_to_csv(file_path, data_list):
    """Write a list of dictionaries to a CSV file"""
    if not data_list:
        return
    
    # Get all possible field names
    fieldnames = set()
    for item in data_list:
        fieldnames.update(item.keys())
    
    # Write to CSV file
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=sorted(fieldnames))
        writer.writeheader()
        writer.writerows(data_list)

def generate_kde_data(ratio_data):
    """Generate KDE data for each interaction type"""
    # Create evaluation points
    x = np.linspace(0, 1, 1000)
    
    # Initialize KDE results
    kde_results = []
    for i in range(len(x)):
        row = {'x': x[i]}
        kde_results.append(row)
    
    # Calculate KDE for each interaction type
    for col, ratios in ratio_data.items():
        if not ratios:
            # If there is no data, fill with zeros
            for row in kde_results:
                row[col] = 0
            continue
        
        # Calculate KDE
        try:
            kde = stats.gaussian_kde(ratios)
            kde_values = kde(x)
            
            # Add KDE values to the results
            for i, row in enumerate(kde_results):
                row[col] = kde_values[i]
        except Exception as e:
            print(f"Warning: Unable to calculate KDE for {col}: {str(e)}")
            for row in kde_results:
                row[col] = 0
    
    return kde_results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python interaction_analysis.py <input_file.csv>")
        sys.exit(1)
    
    sys.exit(analyze_interaction_data(sys.argv[1]))    