import os
from typing import Tuple, Dict


def check_gc_content(seq:str, gc_bounds: Tuple[int, int] = (0, 100)) -> bool:
    '''
    Compute GC-content of the DNA seqence and check if it fits into bounds.

    Arguments:
    - seq (str): DNA sequence
    - gc_bounds (tuple or int, default = (0, 100)): GC-content bounds in percent.
    Tuple if contains lower and upper bounds, int if only 
    contains an upper bound.

    Returns: boolean
    '''
    gc_count= 0
    for nucleotide in seq:
        if nucleotide == 'G' or nucleotide == 'C':
            gc_count += 1
    gc_content = (gc_count /len(seq)) * 100

    if isinstance(gc_bounds, int) or  isinstance(gc_bounds, float):
        gc_bounds = (0, gc_bounds)
        
    return gc_bounds[0] <= gc_content <= gc_bounds[1]

def check_length(seq:str, length_bounds: Tuple[int, int] = (0, 2**32))-> bool:
    '''
    Check if length of the input DNA seqence fits into bounds.

    Arguments:
    - seq (str): DNA sequence
    - length_bounds (tuple or int, default = (0, 2**32)): corresponds to its name. 
    Tuple if contains lower and upper bounds, int if only
    contains an upper bound.

    Returns: boolean
    '''
    if isinstance(length_bounds, int) or  isinstance(length_bounds, float):
        length_bounds = (0, length_bounds)
        
    return length_bounds[0] <= len(seq) <= length_bounds[1]

def check_quality(seq:str, quality_threshold:int = 0)-> bool:
    '''
    Check if average phred quality score of sequence exceeds the threshold.

    Arguments:
    - seq (str): sequence of quality-scores (ASCII encoded)
    - quality_threshold (int, default = 0): corresponds to its name

    Returns: boolean
    '''
    sum_quality = 0
    for symbol in seq:
        sum_quality += ord(symbol) - 33
    mean_quality = sum_quality/len(seq)
    
    return mean_quality > quality_threshold


def read_fastq_file(input_path: str) -> Dict[str, Tuple[str, str, str]]:
    """
    Read FASTQ-file

    Arguments:
    - input_path (str): the path to the FASTQ-file 

    Returns:
    - Dict[str, Tuple[str, str]]: dictionary with sequences, 
    where keys -     sequence identifiers (str), and values - 
    a tuples of two strings: sequence and quality.
    """
    with open(input_path, 'r') as fastq_file:
        names = []
        seqs = []
        comments = []
        qualities = []
        fastqs = dict() 
        for line in fastq_file:
            if line.startswith('@'):
                name = line.strip('\n')
                seq = fastq_file.readline().strip('\n')
                seqs.append(seq)
                comment = fastq_file.readline().strip('\n')
                comments.append(comment)
                quality = fastq_file.readline().strip('\n')
                qualities.append(quality)
                fastqs[name] =  (seqs, comments, qualities)

    return fastqs


def write_fastq_file(filtered_seqs: Dict[str, Tuple[str, str, str]], output_filename: str):
    """
    Write results of FASTQ-filtration to the file

    Arguments:
    - filtered_seqs (Dict[str, Tuple[str, str]]): а dictionary 
    with filtered FASTQ-sequencies
    - output_filename (str): name of the output file
    """
    if not os.path.isdir("fastq_filtrator_results"):
        os.mkdir("fastq_filtrator_results")

    with open(f'fastq_filtrator_results/{output_filename}.fastq', 'w') as output_file:
        for name, values in filtered_seqs.items():
            output_file.write(name + '\n')
            output_file.write(values[0] + '\n')
            output_file.write(values[1] + '\n')
            output_file.write(values[2] + '\n')
            