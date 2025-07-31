import argparse
import pandas as pd
import numpy as np
from scipy import stats
import os

def perform_kde(column_data, bandwidth, ignore_zero):
    """Perform KDE analysis on the given column data"""
    # Convert the column data to numeric and coerce errors to NaN
    valid_data = pd.to_numeric(column_data, errors='coerce').dropna()
    if ignore_zero:
        valid_data = valid_data[valid_data != 0]

    if len(valid_data) < 2:
        print(f"Warning: Insufficient valid samples in column data, KDE analysis cannot be performed")
        return None, None

    # Create a KDE object
    kde = stats.gaussian_kde(valid_data, bw_method=bandwidth)

    # Generate evaluation points
    min_val, max_val = valid_data.min(), valid_data.max()
    x = np.linspace(min_val - 0.5 * (max_val - min_val),
                    max_val + 0.5 * (max_val - min_val),
                    1000)

    # Calculate KDE values
    y = kde(x)

    return x, y

def main():
    # Set command-line arguments
    parser = argparse.ArgumentParser(description='Perform KDE analysis on specified columns of a CSV file')
    parser.add_argument('-i', '--input', required=True, help='Path to the input CSV file')
    parser.add_argument('-o', '--output', required=True, help='Path to the output CSV file')
    parser.add_argument('-b', '--bandwidth', type=float, help='KDE bandwidth parameter')
    parser.add_argument('-c', '--columns', nargs='+', type=int, required=True, help='Indices of columns to analyze (0-based)')
    parser.add_argument('--ignore_zero', action='store_true', default=False, help='Ignore zero values in the data')
    parser.add_argument('--delimiter', default=',', help='Delimiter used in the input CSV file. Default is comma (,). Use \\t for tab.')
    parser.add_argument('--read_title', action='store_true', default=True, help='If true, the first row of the CSV is also data, not a header.')

    args = parser.parse_args()

    try:
        # Read the CSV file
        if args.read_title:
            df = pd.read_csv(args.input, delimiter=args.delimiter, header=None, low_memory=False)
        else:
            df = pd.read_csv(args.input, delimiter=args.delimiter, low_memory=False)

        # Convert all columns to numeric and coerce errors to NaN
        df = df.apply(pd.to_numeric, errors='coerce')

        # Check if the specified columns exist
        for col_index in args.columns:
            if col_index >= df.shape[1]:
                raise ValueError(f"Column index {col_index} is out of range. The CSV file has {df.shape[1]} columns.")

        output_data = {}

        for col_index in args.columns:
            if args.read_title:
                col_name = f"Column_{col_index}"
            else:
                col_name = df.columns[col_index]
            col_data = df[col_index] if args.read_title else df[col_name]

            # Perform KDE analysis
            print(f"Performing KDE analysis on column '{col_name}'...")
            col_x, col_y = perform_kde(col_data, args.bandwidth, args.ignore_zero)

            # Calculate mean and standard deviation
            valid_data = pd.to_numeric(col_data, errors='coerce').dropna()
            if args.ignore_zero:
                valid_data = valid_data[valid_data != 0]
            if len(valid_data) > 0:
                mean_val = valid_data.mean()
                std_val = valid_data.std()
                print(f"Column '{col_name}': Mean = {mean_val:.4f}, Standard Deviation = {std_val:.4f}")
            else:
                print(f"Column '{col_name}': No valid data to calculate mean and standard deviation.")

            # Prepare output data
            if col_x is not None and col_y is not None:
                output_data[f'{col_name}_x'] = col_x
                output_data[f'{col_name}_y'] = col_y

        if not output_data:
            raise ValueError("No valid KDE analysis results to save")

        # Create the output directory (if it doesn't exist)
        output_dir = os.path.dirname(args.output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Save the results to a CSV file
        output_df = pd.DataFrame(output_data)
        output_df.to_csv(args.output, index=False)

        print(f"KDE analysis results have been saved to: {args.output}")

    except FileNotFoundError:
        print(f"Error: Input file '{args.input}' not found")
    except pd.errors.ParserError:
        print(f"Error: Unable to parse CSV file '{args.input}', the format may be incorrect")
    except Exception as e:
        print(f"Error: An unknown error occurred - {str(e)}")

if __name__ == "__main__":
    main()
