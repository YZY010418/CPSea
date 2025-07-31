import json
import csv
import os

def process_templates(template_path, csv_path, output_dir="output"):
    """
    Process template JSON and CSV files to generate new JSON files
    
    Parameters:
        template_path: Path to the template JSON file
        csv_path: Path to the CSV file
        output_dir: Output directory, default is "output"
    """
    os.makedirs(output_dir, exist_ok=True)
    
    with open(template_path, 'r') as f:
        template = json.load(f)
    
    with open(csv_path, 'r') as f:
        reader = csv.reader(f)
        
        for row in reader:
            if len(row) < 2:
                print(f"Warning: Skipping incomplete row {row}")
                continue
            
            base_name = row[0].strip()
            sequence = row[1].strip()
            seq_length = len(sequence)
            
            new_json = template.copy()
            new_json["name"] = base_name
            
            if new_json.get("sequences") and len(new_json["sequences"]) > 0:
                new_json["sequences"][0]["protein"]["sequence"] = sequence
            
            if new_json.get("bondedAtomPairs") and len(new_json["bondedAtomPairs"]) > 0:
                if len(new_json["bondedAtomPairs"][0]) > 1:
                    new_json["bondedAtomPairs"][0][1][1] = seq_length
            
            output_path = os.path.join(output_dir, f"{base_name}.json")
            with open(output_path, 'w') as out_f:
                json.dump(new_json, out_f, indent=2)
            
            print(f"Generated file: {output_path}")

if __name__ == "__main__":
    template_file = "template_headtail.json"
    csv_file = "CPSea_isopep_100.csv"
    
    process_templates(template_file, csv_file)
    print("Processing completed!")