#!/bin/bash

# Set the directory containing the files to process
data_dir="./data/fuso_genomes"

# Loop through all .fasta files in the directory
for fasta_file in "$data_dir"/*.fasta
do
  echo "PNA $name"  
  # Strip the file extension and path to get the name
  name=$(basename "$fasta_file" .fasta)

  # Run the command with the extracted name
  sh mason.sh -f "$data_dir/$name.fasta" -g "$data_dir/$name.gff3" -m 1 -i "$name" -p "./data/PNA_valentina.fasta"
done
