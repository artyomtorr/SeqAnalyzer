# SeqAnalyzer

**SeqAnalyzer** provides a set of Python scripts and utilities for processing and analyzing biological sequences. Whether you're filtering FASTQ files, running sequence predictions, or examining sequence characteristics, SeqAnalyzer offers the functionality you need for your bioinformatics tasks.
### Main scripts:
- **`seqanalyzer.py`**: A module containing utilities for filtering FASTQ sequences and representing biological sequences such as nucleic acids and amino acids.
- **`bio_files_processor.py`**: A module for processing sequence files, including functions to convert multi-line FASTA files to single-line format and to extract genes from GenBank files.
- **`custom_random_forest.py`**: A custom implementation of a random forest algorithm with multiprocessing support.
### Additional Resources:
- **`Showcases.ipynb`**: Jupyter notebook showcasing example usage of SeqAnalyzer tools.
- **`test_seqanalyzer.py`**: Script containing unit tests for the functionality provided in `seqanalyzer.py`.

## Contents

### bio_files_processor.py
- `convert_multiline_fasta_to_oneline`: Converts a multi-line FASTA file to a single-line format.
- `select_genes_from_gbk_to_fasta`: Selects genes from a GenBank file that are adjacent to specified genes of interest.
- `OpenFasta`: Custom context manager for reading FASTA files.
### seqanalyzer.py
- `filter_fastq`: Filters sequences in a FASTQ file based on specified criteria.
- `telegram_logger`: Decorator function for logging function execution and sending Telegram messages.
- `run_genscan`: Python API  for GENSCAN prediction the locations and exon-intron structures of genes.
- `RNASequence`, `DNASequence`, `AminoAcidSequence`: Classes for representing nucleic acid and amino acid sequences.

For detailed usage and examples, refer to the documentation in each script.

## Dependencies
- Python 3.x
- Required Python packages are listed in `requirements.txt`.

## Installation

To get the tool **SeqAnalyzer** clone the git repository::
```bash
git clone https://git@github.com:artyomtorr/SeqAnalyzer.git && cd SeqAnalyzer
```

## Contacts
Please use contacts below to reach out with any comments, concerns, or discussions regarding **SeqAnalyzer.** <br>
- *Artyom Toropov* ([Git-Hub](https://github.com/artyomtorr/), [e-mail](toropov.01@bk.ry))
