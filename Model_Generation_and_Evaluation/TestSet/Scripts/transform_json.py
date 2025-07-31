import json
import os
import argparse

def convert_format(motif):
    """Convert motif string to the required format"""
    positions = []
    parts = motif.split('-')
    for part in parts:
        if part.startswith('R'):
            try:
                position = int(part[1:])
                positions.append(["R", [position, " "]])
            except ValueError:
                print(f"Invalid position in motif: {part}")
    return positions

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Convert and format receptor epitope data to JSON files.')
    parser.add_argument('-i', '--input', required=True, help='Path to the input receptor_epitopes.json file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output directory')
    args = parser.parse_args()

    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    # Read input JSON file
    try:
        with open(args.input, 'r') as file:
            receptor_epitopes = json.load(file)
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # Process each entry
    for key, value in receptor_epitopes.items():
        converted_motif = convert_format(value["motif"])
        
        # Ensure data format meets requirements
        formatted_data = []
        for item in converted_motif:
            if isinstance(item, list) and len(item) == 2 and isinstance(item[0], str) and isinstance(item[1], list):
                formatted_data.append(item)
        
        # Build output file path
        output_filename = f"{key}_complex_pocket.json"
        output_path = os.path.join(args.output, output_filename)
        
        # Save formatted data to JSON file
        try:
            with open(output_path, 'w') as file:
                json.dump(formatted_data, file, separators=(',', ':'))
            print(f"Processed {key}. Saved to {output_path}")
        except Exception as e:
            print(f"Error saving file for {key}: {e}")

    print(f"All files have been saved to the '{args.output}' directory.")

if __name__ == "__main__":
    main()
