# SeqAnalyzer
**SeqAnalyzer** is a new tool for bioinformatic analysis, combining 3 main functions:

- ***Nucleotide sequence analysis:*** get reversed, complementary, transcribed and other versions of nucleotide sequences.
- ***Protein sequence analysis:*** get the key characteristics of the protein sequences: molecular weight, hydrophobicity, length and more.
- ***FASTQ file filtration:***  filter FASTQ-sequences depending on quality, sequence length and GC-content.

  
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
- `get_protein_gene`- returns possible variants of DNAs for a given protein sequence
- `count_amino_acids` - calculates the number of each aminoacid in protein sequence
  
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
  - `transcribe` - transcribes DNA sequence to RNA sequence
  - `reverse` - reverses nucleotide sequence
  - `complement` - returns complementary nucleotide sequence
  - `reverse_complement` - returns reversed complementary nucleotide sequence

#### Examples
```{python}
dna_rna_tools('ATG', 'transcribe') # 'AUG'
dna_rna_tools('ATG', 'reverse') # 'GTA'
dna_rna_tools('AtG', 'complement') # 'TaC'
dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```
### FASTQ filtration




