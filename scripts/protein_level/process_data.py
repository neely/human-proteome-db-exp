# modeified to take in sage lfq protein export and make ready for rank sum

import pandas as pd
import sys
import numpy as np

def process_data(input_file, output_file):
    # Read the input TSV file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t')

    # Select columns with headers ending in "mzML.gz"
    selected_columns = df.filter(like='mzML.gz')

    # Find rows with at least one non-zero value in selected columns
    non_zero_rows = df[df[selected_columns.columns].astype(bool).any(axis=1)]

    # Keep the first column (variable identifiers) corresponding to non-zero rows
    first_column = non_zero_rows.iloc[:, 0]

    # Combine the first column and selected data
    selected_data = pd.concat([first_column, non_zero_rows[selected_columns.columns]], axis=1)

    # Transpose the DataFrame and set the first column as new column headers
    transposed_df = selected_data.transpose()
    transposed_df.columns = transposed_df.iloc[0]

    # Remove the first row (duplicate of new column headers)
    transposed_df = transposed_df.iloc[1:]

    # Take the inverse log2 of all values, but comment out if you want to use log2 transformed vals
    transposed_df = np.power(2, transposed_df.astype(float))

    # Add a "Groups" column based on the first column names
    transposed_df['Groups'] = transposed_df.index.str.startswith('B.naive').astype(int)

    # Rename the index to 0 for 'B.naive' and 1 for 'T4.naive'
    transposed_df.index = transposed_df['Groups']

    # Drop the 'Groups' column
    transposed_df = transposed_df.drop(columns=['Groups'])

    # Export the DataFrame to a CSV file
    transposed_df.to_csv(output_file)

    # Export the number of columns to a text file
    num_columns = len(transposed_df.columns)
    with open('num_columns.txt', 'w') as num_columns_file:
        num_columns_file.write(str(num_columns))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python process_data.py input_file.tsv output_file.csv")
    else:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        process_data(input_file, output_file)
