# De-Bruijn-graph-generator
A Python tool to construct De Bruijn graphs from FASTQ files for genomic sequence assembly and short-read analysis.

## Overview
This code builds a **De Bruijn graph** with a fastq file provided by the user.
A De Bruijn graph is a directed graph used in bioinformatics to represent overlaps between sequences of symbols and in this case genomic reads.
By breaking short sequencing reads into overlapping k-mers, this method helps reconstructs longer genomic sequences and chain together the correct nucleotide path.

### Features
- Parses and processes genomic reads from a user-provided FASTQ file.
- Constructs a directed De Bruijn graph using a configurable $k$-mer size.
- Helps visualize and analyze sequence overlaps for genome assembly.

### Pipeline Overview
The 5 main steps of this pipeline are:
1. **FASTQ Parsing:** Reading and extracting raw nucleotide sequences from the input file.
2. **$k$-mer Extraction:** Splitting sequencing reads into overlapping subsequences of length $k$.
3. **Graph Construction:** Building the directed De Bruijn graph from $(k-1)$-mer prefixes and suffixes.
4. **Eulerian Path Traversal:** Finding a path through the graph to reconstruct the contiguous sequence (contig).
5. **Graph Export:** Exporting the network structure into a standard `.dot` file.

> **Tip for Visualization:** You can render and explore the generated `.dot` file directly in VS Code using the **Graphviz Interactive Preview** extension, or paste its contents into an online viewer like [Edotor](https://edotor.net/).
