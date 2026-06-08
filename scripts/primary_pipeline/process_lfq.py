import os
import pandas as pd

# Function to read and filter lfq.tsv files
def process_lfq_file(file_path):
    df = pd.read_csv(file_path, delimiter='\t')
    filtered_df = df[df['q_value'] < 0.01]  # Update column name to 'q_value'
    return set(filtered_df['peptide'])

# Function to combine and remove duplicates from lists
def combine_and_remove_duplicates(lists):
    combined_list = []
    for peptide_list in lists:
        combined_list.extend(peptide_list)
    return list(set(combined_list))

# Function to write CSV and output.txt files
def write_output_files(unique_peptides):
    with open('output.txt', 'w') as output_file:
        output_file.write(f'Number of Unique Peptide IDs: {len(unique_peptides)}')

    df = pd.DataFrame({'Peptide ID': unique_peptides})
    df.to_csv('output.csv', index=False)

# Main script
folder_path = '.'  # Current directory

peptide_lists = []

for folder_name in os.listdir(folder_path):
    file_path = os.path.join(folder_path, folder_name, 'lfq.tsv')
    if os.path.isfile(file_path):
        peptide_list = process_lfq_file(file_path)
        peptide_lists.append(peptide_list)

unique_peptides = combine_and_remove_duplicates(peptide_lists)
write_output_files(unique_peptides)
