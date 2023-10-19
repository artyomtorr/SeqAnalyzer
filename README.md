# SeqAnalyzer
**SeqAnalyzer** is a new tool for bioinformatic analysis, combining 3 main functions:

- ***Nucleotide sequence analysis:*** get reversed, complementary, transcribed and other versions of nucleotide sequences.
- ***Protein sequence analysis:*** get the key characteristics of the protein sequences: molecular weight, hydrophobicity, length and more.
- ***FASTQ file filtration:***  filter FASTQ-sequences based on quality, sequence length and GC-content.
- ***Processing of GBK and FASTA files***: convert FASTA file from multiline format to oneline or select genes from GBK file that are adjacent to the genes of interest.

## Installation

To get the tool **SeqAnalyzer** clone the git repository::
```bash
git clone https://git@github.com:artyomtorr/SeqAnalyzer.git && cd SeqAnalyzer
```
## Usage
### Nucleotide and protein sequences analsysis
To start working with nucleotide or protein sequences, run appropriate **command** (see below), which takes arbitrary number of arguments with protein or nucleotide sequencies (*str*) and the name of the procedure to be performed (always the last argument, *str*):
```{python}
run_dna_rna_tools(*sequences: str, procedure: str) #Nucleotide sequences analsysis
run_protein_tools(*sequences: str, procedure: str) #Protein sequences analsysis
```
**NOTE:**  The procedure `check_mutations` for protein sequences has a fixed number of string arguments: one RNA sequence, one protein sequence and the name of procedure itself.

The **output** is *dictionary* with protein sequences as keys and the results of the procedure as values. If more than one sequence is submitted the output is list of dictionaries of the same format.

### Protein procedures

- `compute_molecular_weight` — computes molecular weight of protein sequence in g/mol
- `compute_length` — computes the number of amino acids in protein sequence
- `compute_hydrophobicity` — computes the percentage of gydrophobic aminoacids in protein sequence
- `check_mutations` — checks missense mutations in the protein sequence after translation
- `get_protein_gene`— returns possible variants of DNAs for a given protein sequence
- `count_amino_acids` — calculates the number of each aminoacid in protein sequence
  
#### Examples
```{python}
run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_length')
#[{'MAEGEITNLP': 10}, {'tGQYLAMDTSgLLYGSQT': 18}]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_molecular_weight')
#[{'MAEGEITNLP': 1055.496}, {'tGQYLAMDTSgLLYGSQT': 1886.872}]

run_protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_hydrophobicity')
#[{'MAEGEITNLP': 50.0}, {'tGQYLAMDTSgLLYGSQT': 27.778}]

run_protein_tools('AUGGAUCAUcAAUAA', 'MDKL*', 'check_mutations')
#{'MDKL*': 'Mutations: K3, L4.'}

run_protein_tools('MAEGLP', 'get_protein_gene')
#{'MAEGLP': 'ATG GCT/GCC/GCA/GCG GAA/GAG GGT/GGC/GGA/GGG TTA/TTG/CTT/CTC/CTA/CTG CCT/CCC/CCA/CCG'}

run_protein_tools('MAEGLP', 'LYGSQT','stats_amino_acids')
#[{'MAEGLP': {'M': 1, 'A': 1, 'E': 1, 'G': 1, 'L': 1, 'P': 1}},
#{'LYGSQT': {'L': 1, 'Y': 1, 'G': 1, 'S': 1, 'Q': 1, 'T': 1}}]
```
### Nucleotide sequences procedures
  - `transcribe` — transcribes DNA sequence to RNA sequence
  - `reverse` — reverses nucleotide sequence
  - `complement` — returns complementary nucleotide sequence
  - `reverse_complement` — returns reversed complementary nucleotide sequence

#### Examples
```{python}
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```
### FASTQ filtration
If you want to filter FASTQ-sequences based on quality, sequence length and GC-content, run this **command**:
```{python}
run_fastq_filter(input_path: str, output_filename: str, gc_bounds = (0, 100), length_bounds = (0, 2**32), quality_threshold = 0) 
```
In exapmle you can see few arguments:
- `input_path`(str): the path to the file with FASTQ-sequences
- `output_filename` (str): the name of the output file with filtered FASTQ-sequences
- `gc_bounds` (tuple or int, default = (0, 100)): GC-content interval (percentage) for filtering. Tuple if contains lower and upper bounds, int if only contains an upper bound.
- `length_bounds` (tuple or int, default = (0, 2**32)): length interval for filtering. Tuple if contains lower and upper bounds, int if only contains an upper bound.
- `quality_threshold` (int, default = 0): threshold value of average read quality for filtering.

The **output** FASTQ-file consists sequences that meet all the requirements. All files with filtration results are stored in the the */fastq_filtrator_results*
directory. The name of the input file is used as default output file name.

### Processing of GBK and FASTA files
**bio_files_processor.py** - tool for for processing of GBK and FASTA files. <br>
If you want to *convert FASTA file from multiline format to oneline*, run this **command**:
```{python}
convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None)
```
In exapmle you can see few arguments:
- `input_fasta` (str): name of the input multiline FASTA file.
- `output_fasta` (str): name of the output oneline FASTA file. The name of the input file is used as default output file name. <br>

As a **result** of execution, the function creates a FASTA file with same sequences but in oneline format.

If you want to *select genes from GBK file that are adjacent to the genes of interest*, run this **command**:
```{python}
select_genes_from_gbk_to_fasta(input_gbk: str, genes: list, n_before: int = 1, n_after: int = 1, output_fasta: str = None):
```
In exapmle you can see few arguments:
- `input_gbk` (str): the path to the input GBK file.
- `genes` (list): a list of names to genes of interest.
- `n_before` (int, default = 1): number of genes to extract before the gene of interest.
- `n_after` (int, default = 1): number of genes to extract after the gene of interest.
- `output_fasta` (str, optional): name of the output FASTA file. The name of the input file is used as default output file name.

As a **result** of execution, the function creates a FASTA file with genes that are adjacent to the genes of interest and their translation sequences.
### Troubleshooting
- The tool works **only** with protein and RNA/DNA sequences. If any of the entered sequences contain inappropriate characters or cannot exist, the program will display an error. Sequences can contain characters of any case.

```python
run_protein_tools('PROTEIN', 'compute_molecular_weight') #ValueError: Invalid protein sequence
run_protein_tools('AUGGAU_AUcAAUAA', 'MDKL*', 'check_mutations') #ValueError: Invalid RNA sequence
```
- For the protein procedure `check_mutations` there are extra requirements for RNA and protein sequences: mRNA sequences must contain **start-codon** and **one of the stop-codons**, protein sequnces must start with **"M"** and ends with **"*"** (stop-codon). 
```python
run_protein_tools("AUGGUAGGGAAAUUUUGA", "MGGKF", 'check_mutations') #ValueError: Stop (*) is absent
run_protein_tools("AUGGUAGGGAAAUUUUGA", "GGKF*", 'check_mutations') #ValueError: Start (M) is absent
```
### Contacts
Please use contacts below to reach out with any comments, concerns, or discussions regarding **SeqAnalyzer.** <br>
- *Artyom Toropov* ([Git-Hub](https://github.com/artyomtorr/), [e-mail](toropov.01@bk.ry))


