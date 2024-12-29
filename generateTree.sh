#!/bin/bash

# Check if at least one file was provided as an argument
if [ "$#" -lt 1 ]; then
    echo "Error: No files provided. Please provide one or more fasta files as arguments."
    exit 1
fi

# Loop through each file passed as a command-line argument
for file in "$@"
do
   echo "Running toxins for $file"
   
   # Extract toxins
   awk -v RS=">" -v FS="\n" -v ORS="" '$1 ~ /\.toxin$/ {print ">"$0}' "$file" > "$file.new.toxin.fasta"
   muscle -align "$file.new.toxin.fasta" -output "$file.new.toxin.afa"
   FastTree < "$file.new.toxin.afa" > "$file.new.toxin.nwk"

   # Repeat for antitoxins
   echo "Running antitoxins for $file"
   awk -v RS=">" -v FS="\n" -v ORS="" '$1 ~ /\.antitoxin$/ {print ">"$0}' "$file" > "$file.new.antitoxin.fasta"
   muscle -align "$file.new.antitoxin.fasta" -output "$file.new.antitoxin.afa"
   FastTree < "$file.new.antitoxin.afa" > "$file.new.antitoxin.nwk"
   
   # Just for readability
   echo ""
done

