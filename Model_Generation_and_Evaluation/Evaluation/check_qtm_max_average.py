import sys
import csv
import argparse

def process_tsv_file(file_path, ignore_R, ignore_same_id, overlap_threshold, column):
    max_values = {}  # key: first column value, value: (max value, entire row data)
    
    try:
        with open(file_path, 'r', encoding='utf-8', newline='') as file:
            reader = csv.reader(file, delimiter='\t')
            
            for row_num, row in enumerate(reader, 1):
                if row_num % 1000000 == 0:
                    print(f"Processed {row_num:,} rows...", file=sys.stderr)
                
                if len(row) < column:
                    print(f"Warning: Row {row_num} has fewer than {column} columns, skipped", file=sys.stderr)
                    continue
                
                key = row[0]
                second_col = row[1] if len(row) > 1 else ""
                
                if ignore_same_id:
                    key_prefix = key.split('_', 1)[0] if '_' in key else key
                    key_prefix_overlap = key_prefix[:overlap_threshold]
                    
                    second_col_prefix = second_col.split('-', 1)[0] if '-' in second_col else second_col
                    second_col_prefix_overlap = second_col_prefix[:overlap_threshold]
                    
                    if key_prefix_overlap == second_col_prefix_overlap:
                        continue
                
                if ignore_R:
                    last_underscore_idx = key.rfind('_')
                    if last_underscore_idx != -1 and last_underscore_idx < len(key) - 1:
                        suffix = key[last_underscore_idx + 1:]
                        if suffix == 'R':
                            continue
                
                try:
                    value_str = row[column - 1].replace('^@', '').strip()
                    value = float(value_str)
                except ValueError:
                    print(f"Warning: Row {row_num} column {column} '{row[column - 1]}' is not a valid number, skipped", file=sys.stderr)
                    continue
                
                if key not in max_values or value > max_values[key][0]:
                    max_values[key] = (value, row)
        
        all_max_values = [item[0] for item in max_values.values()]
        if not all_max_values:
            print("Error: No valid numerical data found", file=sys.stderr)
            return
        
        average = sum(all_max_values) / len(all_max_values)
        
        print("----------------------------------------")
        print("\n----------------------------------------")
        print(f"Average of maximum values: {average:.6f}")
        print(f"Number of groups: {len(all_max_values)}")
        
    except FileNotFoundError:
        print(f"Error: File '{file_path}' does not exist", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error processing file: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process TSV file, group by first column and find maximum value of specified column in each group')
    parser.add_argument('file_path', help='Path to TSV file')
    parser.add_argument('--ignore_R', action='store_true', default=False, 
                       help='Ignore rows where the part after the last underscore in first column is R (default: do not ignore)')
    parser.add_argument('--ignore_same_id', action='store_true', default=False,
                       help='Ignore rows where specified number of prefix characters are the same between first column (split by first _) and second column (split by first -) (default: do not ignore)')
    parser.add_argument('--overlap_threshold', type=int, default=4,
                       help='When --ignore_same_id is true, compare first N characters (default: 4)')
    parser.add_argument('-c', '--column', type=int, choices=[3, 4], default=3,
                       help='Specify whether to statistics on column 3 or 4 (default: 3)')
    
    args = parser.parse_args()
    
    if args.overlap_threshold <= 0:
        print("Error: --overlap_threshold must be a positive integer", file=sys.stderr)
        sys.exit(1)
    
    print(f"Starting to process file: {args.file_path}")
    print(f"Statistics column: Column {args.column}")
    print(f"Ignore entries ending with _R: {'Yes' if args.ignore_R else 'No'}")
    print(f"Ignore entries with same prefix: {'Yes' if args.ignore_same_id else 'No'}")
    if args.ignore_same_id:
        print(f"Compare first {args.overlap_threshold} characters")
    
    process_tsv_file(args.file_path, args.ignore_R, args.ignore_same_id, args.overlap_threshold, args.column)