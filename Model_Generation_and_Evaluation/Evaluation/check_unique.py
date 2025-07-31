import argparse

def extract_unique_first_column(input_path, output_path, encoding='utf-8'):
    """
    Extract non-redundant content from the first column of a tab-separated file
    
    Parameters:
        input_path: Path to the input TSV file
        output_path: Path to the output file (one unique value per line)
        encoding: File encoding format, default is utf-8
    """
    unique_values = set()
    
    try:
        with open(input_path, 'r', encoding=encoding) as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.rstrip('\n')
                first_col = line.split('\t', 1)[0]
                unique_values.add(first_col)
                
                if line_num % 100000 == 0:
                    print(f"Processed {line_num} lines, current unique values count: {len(unique_values)}")
        
        with open(output_path, 'w', encoding=encoding) as outfile:
            for value in unique_values:
                outfile.write(f"{value}\n")
        
        print(f"Processing completed! Total lines processed: {line_num}, final unique values count: {len(unique_values)}")
        print(f"Results saved to: {output_path}")
    
    except Exception as e:
        print(f"Error occurred: {str(e)}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Extract unique values from the first column of a TSV file')
    parser.add_argument('-i', '--input', required=True, help='Input TSV file path')
    parser.add_argument('-o', '--output', required=True, help='Output file path')
    args = parser.parse_args()
    
    extract_unique_first_column(args.input, args.output)