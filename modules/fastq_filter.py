from typing import Tuple, Dict
import os


def check_gc_content(seq:str, gc_bounds = (0, 100)) -> bool:
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

    if isinstance(gc_bounds, tuple):
            lower_bound = gc_bounds[0]
            upper_bound = gc_bounds[1]
            return gc_content >= lower_bound and gc_content <= upper_bound
    else:
        return gc_content < gc_bounds

def check_length(seq:str, length_bounds = (0, 2**32))-> bool:
    '''
    Check if length of the input DNA seqence fits into bounds.

    Arguments:
    - seq (str): DNA sequence
    - length_bounds (tuple or int, default = (0, 2**32)): corresponds to its name. 
    Tuple if contains lower and upper bounds, int if only
    contains an upper bound.

    Returns: boolean
    '''
    if isinstance(length_bounds, tuple):
            lower_bound = length_bounds[0]
            upper_bound = length_bounds[1]
            return len(seq) >= lower_bound and len(seq) <= upper_bound
    else:
        return len(seq) < length_bounds


def check_quality(seq:str, quality_threshold = 0)-> bool:
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
    quality = sum_quality/len(seq)
    
    return quality > quality_threshold


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
        for line in fastq_file:
            if line.startswith('@SRX'):
                name = line.strip('\n')
                names.append(name)
                seq = fastq_file.readline().strip('\n')
                seqs.append(seq)
                comment = fastq_file.readline().strip('\n')
                comments.append(comment)
                quality = fastq_file.readline().strip('\n')
                qualities.append(quality)
        keys = names
        values = list(zip(seqs, comments, qualities))
        fastq_dict = dict(zip(keys, values))
    return fastq_dict


def write_fastq_file(filtered_seqs: Dict[str, Tuple[str, str, str]], output_filename):
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
            