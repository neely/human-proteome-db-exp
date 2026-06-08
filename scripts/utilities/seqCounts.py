import os
import csv

# List all .fasta files in the current directory
fasta_files = [file for file in os.listdir() if file.endswith('.fasta')]

# Create a CSV file for output
with open('output.csv', 'w', newline='') as csvfile:
    csvwriter = csv.writer(csvfile)

    # Write header
    csvwriter.writerow(['Filename', 'Count'])

    # Iterate through each .fasta file
    for file in fasta_files:
        with open(file, 'r') as fasta:
            # Count occurrences of '>'
            count = sum(1 for line in fasta if line.startswith('>'))

            # Write filename and count to CSV
            csvwriter.writerow([file, count])
