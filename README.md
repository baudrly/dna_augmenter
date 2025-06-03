# DNA Sequence Augmenter

A Rust program for generating augmented DNA sequences from FASTA files using various strategies.

## Features

- Multiple augmentation strategies:
  - Sliding windows
  - Sequence chunks
  - Centered windows
  - Truncations (beginning, end, random)
  - Synonymous codon variants
- Reverse complement generation (for coding sequences only)
- Custom genetic code support
- Parallel processing
- FASTA format input/output

## Installation

1. Clone this repository
2. Build with `cargo build --release`

## Usage

```bash
dna_augmenter -i input.fa -o output.fa [OPTIONS] --strategy <STRATEGY>...
```

### Required Arguments

```
-i, --input-path <FILE>: Input FASTA file

-o, --output-path <FILE>: Output FASTA file

--strategy <STRATEGY>: One or more augmentation strategies (see below)
```

### Optional Arguments

```
--revcomp: Also generate reverse complements of all augmented sequences

--genetic-code-file <FILE>: Custom genetic code definition file

-t, --threads <NUM>: Number of threads to use (default: auto-detect)
```

### Strategy Format

Strategies are specified in the format: `type:count[:param1:param2]`. 

* **Sliding windows**: `sliding:count:window_size:step_size`

*Example*: `sliding:5:300:100` - 5 windows of size 300 with step 100

* **Chunks**: `chunks:count:num_chunks`

*Example*: `chunks:2:3` - 2 sets of 3 equal chunks

* **Centered window**: `centered:count:window_size`

*Example*: `centered:3:500` - 3 centered windows of 500bp

* **Beginning truncation**: `begin_trunc:count:min_keep_len:max_keep_len`

*Example*: `begin_trunc:4:800:1200` - 4 truncations keeping 800-1200bp from start

* **End truncation**: `end_trunc:count:min_keep_len:max_keep_len`

*Example*: `end_trunc:4:800:1200` - 4 truncations keeping 800-1200bp from end

* **Random truncation**: `random_trunc:count:min_output_len:max_output_len`

*Example*: `random_trunc:5:500:1000` - 5 random subsequences of 500-1000bp

* **Synonymous variants**: `synonymous:count`

*Example*: `synonymous:10` - 10 synonymous codon variants

### Genetic Code File Format

Each line should contain:

```
CODON AMINO_ACID [S]
```

Where:

`CODON`: 3 DNA bases (e.g., "ATG")

`AMINO_ACID`: Single letter code (e.g., "M")

Optional "S" marks start codons

## Examples

```bash
dna_augmenter -i genes.fa -o augmented.fa \
    --strategy sliding:5:300:100 \
    --strategy chunks:2:3 \
    --strategy synonymous:10 \
    --revcomp \
    --threads 8
```

This will:

* Create 5 sliding windows (300bp, step 100bp) for each sequence
* Create 2 sets of 3 chunks for each sequence
* Create 10 synonymous variants for each sequence
* Generate reverse complements for all augmented sequences
* Use 8 threads for processing

## Output

Output sequences will have headers in the format:

```
>original_header_augN_strategy
```

Where:

* N is the augmentation index
* strategy describes the augmentation method

Reverse complements will have _rc appended to the header.