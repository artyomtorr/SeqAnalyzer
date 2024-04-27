import os
from dataclasses import dataclass
from typing import Dict, List, Tuple


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Convert FASTA file from multiline format to oneline.

    Arguments:
    - input_fasta (str): name of the input multiline FASTA file 
    - output_fasta (str): name of the output oneline FASTA file 
    
    The default output file name is the name of the input file.
    """
    seqs = {}
    with open(input_fasta) as fasta_file:
        name = ''
        seq = []
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if name and seq: # every time a new name starts, add full oneline sequence into the dict
                    seqs[name] = ''.join(seq)
                name = line # initialize new name
                seq = [] # reset sequence to empty list
            else:
                seq.append(line)
        if name and seq: # add last sequence
            seqs[name] = ''.join(seq)
            
    if not output_fasta:
        output_fasta = os.path.splitext(input_fasta)[0] 

    with open(f'{output_fasta}.fasta', 'w') as output_file:
        for name, seq in seqs.items():
            output_file.write(name + '\n')
            output_file.write(seq + '\n')


def get_genes_from_gbk(input_gbk: str) -> Dict[str, Tuple[str,str]]:
    """
    Extract genes and ther translation sequences from GBK file

    Arguments:
    - input_gbk(str): the path to the input GBK file 

    Returns:
    - genes_dict(Dict[str, Tuple[str,str]]): a dict with extracted genes
    and their translation sequences in the format:
    {gene: (order number, translation sequence)}
    """
    genes_dict = {}
    gene_number = 0
    index = 0
    seqs = []
    with open(input_gbk, 'r') as file:
        lines = file.readlines()
        while index < len(lines):
            lines[index] = lines[index].strip()
            if lines[index].startswith('/gene='): 
                lines[index] = lines[index].strip('"\n')
                gene_name = lines[index].replace('/gene="', '') 
                index += 1
                while not lines[index].strip().startswith('/translation='):  
                    index += 1
                while not lines[index].endswith('"\n'): #add multiline translation sequence
                    seqs.append(lines[index].strip())
                    index += 1
                seqs.append(lines[index].strip().rstrip('"')) #add oneline or last line of multiline translation sequence
                full_seq = ''.join(seqs).replace('/translation="', '')
                seqs = []
                genes_dict[gene_name] = gene_number, full_seq
                gene_number +=1
            index += 1
    return genes_dict


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list, 
                                   n_before: int = 1, n_after: int = 1,
                                   output_fasta: str = None) -> None:
    """
    Select genes from GBK file that are adjacent to the genes of interest

    Arguments:
    - input_gbk (str): the path to the input GBK file 
    - genes (list): a list of names to genes of interest
    - n_before (int, default = 1): number of genes to extract before the gene of interest
    - n_after (int, default = 1): number of genes to extract after the gene of interest
    - output_fasta (str, optional): name of the output FASTA file 
    The default output file name is the name of the input file.
    """
    genes_dict = get_genes_from_gbk(input_gbk)
    genes_from_gbk = list(genes_dict.keys()) 
    max_number = len(genes_from_gbk)
    genes_of_interest = []

    for gene in genes:
        gene_number = genes_dict.get(gene)[0]
        if gene_number - n_before > 0:
            genes_of_interest_start = gene_number - n_before 
        else:
            genes_of_interest_start = 0
        if gene_number + n_after < max_number:
            genes_of_interest_end = gene_number + n_after
        else:
            genes_of_interest_end = max_number

        genes_of_interest.extend(genes_from_gbk[genes_of_interest_start:gene_number])
        genes_of_interest.extend(genes_from_gbk[gene_number+1:genes_of_interest_end+1])

    genes_of_interest = list(dict.fromkeys(genes_of_interest))  

    if not output_fasta:  
        output_fasta = os.path.basename(input_gbk)
        output_fasta = os.path.splitext(output_fasta)[0]

    with open(f'{output_fasta}.fasta', mode='w') as output_file:
        for gene in genes_of_interest:
            output_file.write(">" + gene + "\n")
            output_file.write(genes_dict.get(gene)[1] + "\n")


@dataclass
class FastaRecord:
    """Dataclass representing a FASTA record."""
    id: str
    description: str
    seq: str

    def __repr__(self)-> str:
        line_length = 75
        seq_lines = [self.seq[i:i+line_length] for i in range(0, len(self.seq), line_length)]
        output = f">{self.id} {self.description}\n"
        output += '\n'.join(seq_lines)
        return output
    

class OpenFasta:
    """Context manager for reading FASTA files."""
    def __init__(self, file_path: str, mode: str = 'r') -> None:
        self.file_path = file_path
        self.mode = mode
        self.file_handler = None
        self.current = None

    def __enter__(self) -> 'OpenFasta':
        self.file_handler = open(self.file_path, mode=self.mode)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb) -> None:
        if self.file_handler:
            self.file_handler.close()

    def __iter__(self) -> 'OpenFasta':
        return self

    def __next__(self) -> FastaRecord:
        return self.read_record()

    def read_record(self) -> FastaRecord:
        record_id = ""
        description = ""
        sequence = ""

        if self.current:
            line = self.current
        else:
            line = self.file_handler.readline().strip()

        if line == "":
            raise StopIteration

        if line.startswith(">"):
            record_id, description = line[1:].split(maxsplit=1)

        while True:
            line = self.file_handler.readline().strip()
            if not line or line.startswith(">"):
                self.current = line
                break
            sequence += line

        return FastaRecord(record_id, description, sequence)

    def read_records(self) -> List[FastaRecord]:
        records = []
        while True:
            try:
                record = self.read_record()
                records.append(record)
            except StopIteration:
                break
        return records
    