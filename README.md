# SeqAnalyzer
**SeqAnalyzer** is a new tool for bioinformatic analysis, combining 3 main functions:

- ***Nucleotide sequence analysis:*** get reversed, complementary, transcribed and other versions of nucleotide sequences.
- ***Protein sequence analysis:*** get the key characteristics of the protein sequences: molecular weight, hydrophobicity, length and more.
- ***FASTQ file filtration:***  filter FASTQ-sequences based on quality, sequence length and GC-content.

  
## Installation

To get the tool **SeqAnalyzer** clone the git repository::
```bash
git clone https://git@github.com:artyomtorr/SeqAnalyzer.git && cd SeqAnalyzer
```
## Usage
### Nucleotide and protein sequences analsysis
To start working with nucleotide or protein sequences, run appropriate **command** (see below), which takes arbitrary number of arguments with protein or nucleotide sequencies (*str*) and the name of the procedure to be performed (always the last argument, *str*):
```{python}
dna_rna_tools(*sequences: str, procedure: str) #Nucleotide sequences analsysis
protein_tools(*sequences: str, procedure: str) #Protein sequences analsysis
```
**NOTE:**  The procedure `check_mutations` for protein sequences has a fixed number of string arguments: one RNA sequence, one protein sequence and the name of procedure itself.

The **output** is the result of the procedure as *string*, or *dictionary* if one sequence is submitted or *list* of strings/dictionaries if several.

### Protein procedures

- `compute_molecular_weight` — computes molecular weight of protein sequence in g/mol
- `compute_length` — computes the number of amino acids in protein sequence
- `compute_hydrophobicity` — computes the percentage of gydrophobic aminoacids in protein sequence
- `check_mutations` — checks missense mutations in the protein sequence after translation
- `get_protein_gene`— returns possible variants of DNAs for a given protein sequence
- `count_amino_acids` — calculates the number of each aminoacid in protein sequence
  
#### Examples
```{python}
protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_length')
#[{'MAEGEITNLP': 10}, {'tGQYLAMDTSgLLYGSQT': 18}]

protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_molecular_weight')
#[{'MAEGEITNLP': 1055.496}, {'tGQYLAMDTSgLLYGSQT': 1886.872}]

protein_tools('MAEGEITNLP', 'tGQYLAMDTSgLLYGSQT', 'compute_hydrophobicity')
#[{'MAEGEITNLP': 50.0}, {'tGQYLAMDTSgLLYGSQT': 27.778}]

protein_tools('AUGGAUCAUcAAUAA', 'MDKL*', 'check_mutations')
#'Mutations: K3, L4.'

protein_tools('MAEGLP', 'LYGSQT','get_protein_gene')
#['ATG GCT/GCC/GCA/GCG GAA/GAG GGT/GGC/GGA/GGG TTA/TTG/CTT/CTC/CTA/CTG CCT/CCC/CCA/CCG',
#'TTA/TTG/CTT/CTC/CTA/CTG TAT/TAC GGT/GGC/GGA/GGG TCT/TCC/TCA/TCG/AGT/AGC CAA/CAG ACT/ACC/ACA/ACG']

protein_tools('MAEGLP', 'LYGSQT','stats_amino_acids')
#[{'M': 1, 'A': 1, 'E': 1, 'G': 1, 'L': 1, 'P': 1},
#{'L': 1, 'Y': 1, 'G': 1, 'S': 1, 'Q': 1, 'T': 1}]
```
### Nucleotide sequences procedures
  - `transcribe` — transcribes DNA sequence to RNA sequence
  - `reverse` — reverses nucleotide sequence
  - `complement` — returns complementary nucleotide sequence
  - `reverse_complement` — returns reversed complementary nucleotide sequence

#### Examples
```{python}
dna_rna_tools('ATG', 'transcribe') # 'AUG'
dna_rna_tools('ATG', 'reverse') # 'GTA'
dna_rna_tools('AtG', 'complement') # 'TaC'
dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```
### FASTQ filtration
If you want to filter FASTQ-sequences based on quality, sequence length and GC-content, run this **command**:
```{python}
fastq_filter(seqs: dict, c_bounds = (0, 100),  length_bounds = (0, 2**32), quality_threshold = 0) 
```
In exapmle you can see few arguments:
- `seqs`: a dict of FASTQ-sequences in the format: 'id' : ('sequence', 'quality'). Where key - sequence identifier (str), and value is a tuple of two strings: sequence and quality.
- `gc_bounds` (tuple or int, default = (0, 100)): GC-content interval (percentage) for filtering. Tuple if contains lower and upper bounds, int if only contains an upper bound.
- `length_bounds` (tuple or int, default = (0, 2**32)): length interval for filtering. Tuple if contains lower and upper bounds, int if only contains an upper bound.
- `quality_threshold` (int, default = 0): threshold value of average read quality for filtering.

The **output**  is a dict consisting FASTQ-sequences that meet all the requirements.

### Troubleshooting
- The tool works **only** with protein and RNA/DNA sequences. If any of the entered sequences contain inappropriate characters or cannot exist, the program will display an error. Sequences can contain characters of any case.

```python
protein_tools('PROTEIN', 'compute_molecular_weight') #ValueError: Invalid protein sequence
protein_tools('AUGGAU_AUcAAUAA', 'MDKL*', 'check_mutations') #ValueError: Invalid RNA sequence
```
- For the protein procedure `check_mutations` there are extra requirements for RNA and protein sequences: mRNA sequences must contain **start-codon** and **one of the stop-codons**, protein sequnces must start with **"M"** and ends with **"*"** (stop-codon). 
```python
run_protein_tools("AUGGUAGGGAAAUUUUGA", "MGGKF", 'check_mutations') #ValueError: Stop (*) is absent
run_protein_tools("AUGGUAGGGAAAUUUUGA", "GGKF*", 'check_mutations') #ValueError: Start (M) is absent
```
### Contacts
Please use contacts below to reach out with any comments, concerns, or discussions regarding **SeqAnalyzer.** <br>
- *Artyom Toropov* ([Git-Hub](https://github.com/artyomtorr/), [e-mail](toropov.01@bk.ry))


