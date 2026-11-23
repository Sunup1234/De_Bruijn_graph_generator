# De-Bruijn-graph-generator
A Python tool to construct De Bruijn graphs from FASTQ files for genomic sequence assembly and short-read analysis.

## Overview
This code builds a **De Bruijn graph** with a fastq file provided by the user.
A De Bruijn graph is a directed graph used in bioinformatics to represent overlaps between sequences of symbols and in this case genomic reads.
By breaking short sequencing reads into overlapping k-mers, this method helps reconstructs longer genomic sequences and chain together the correct nucleotide path.

### Features
- Parses and processes genomic reads from the fastq file provided by the user.
- Constructs a directed De Bruijn graph based on configurable k-mer sizes.
- Helps to vizualize/Analysing sequence overlaps for genome assembly.
