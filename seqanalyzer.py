# Importing modules
import modules.fastq_filter 
import modules.dna_rna_tools 
import modules.protein_tools 

# Function for FASTQ-sequences filtration
def fastq_filter(seqs: dict, gc_bounds = (0, 100), 
                 length_bounds = (0, 2**32), 
                 quality_threshold = 0) -> dict:
    '''
    Filter FASTQ-sequences based on entered requirements.
    
    Arguments:
        - seqs: a dict of FASTQ-sequences in the format:
    'id' : ('sequence', 'quality'). Where key - sequence identifier 
    (str), and value - a tuple of two strings: sequence and quality.
        - gc_bounds (tuple or int, default = (0, 100)): GC-content interval 
    (percentage) for filtering. Tuple if contains lower and upper 
    bounds, int if only contains an upper bound.
        - length_bounds (tuple or int, default = (0, 2**32)): length interval 
    for filtering. Tuple if contains lower and upper bounds, 
    int if only contains an upper bound.
        - quality_threshold (int, default = 0): threshold value of average 
    read quality for filtering.

    Returns: a dict consisting FASTQ-sequences that meet all the requirements.
    '''
    result = {}
    for id in seqs:
        quality_value = check_quality(seqs[id][1], quality_threshold) 
        length_value = check_length(seqs[id][0], length_bounds)
        gc_content_value = check_gc_content(seqs[id][0], gc_bounds)
        if  quality_value and length_value and gc_content_value == True:
            result[id] = seqs[id]
    return result

# Function for DNA/RNA sequences analysis
def dna_rna_tools(*args: str):
    """
    Function containing methods for DNA or RNA sequences analysis.

    Takes arbitrary number of arguments with DNA or RNA sequencies
    and the name of the procedure to be performed (always the last
    argument). 
    Returns the result of the procedure as string if one sequence
    is submitted or list if several.
    """
    *seqs, procedure = args
    results = []
    functions = {'transcribe': transcribe, 
                 'reverse': reverse,
                 'reverse_complement': reverse_complement,
                 'complement': complement}
    for seq in seqs:
        if is_rna(seq) and is_dna(seq) is not True:
            raise ValueError("Invalid alphabet")
        results.append(functions[procedure](seq))
    if len(results) == 1:
        return results[0]
    else: 
        return results

# Function for protein sequences analysis
def protein_tools(*args: str):
    """
    Function containing methods for protein analysis.

    Takes arbitrary number of arguments with protein sequencies
    and the name of the procedure to be performed (always the last
    argument). Returns the result of the procedure as string
    or dictionary if one sequence is submitted or list if several.

    Note: if procedure 'check_mutations' is used then input must
    contain only three arguments: RNA sequence, protein sequence
    and the name of procedure itself.
    """
    *seqs, procedure = args
    results = []
    functions = {
        "compute_molecular_weight": compute_molecular_weight,
        "compute_length": compute_length,
        "compute_hydrophobicity": compute_hydrophobicity,
        "stats_amino_acids": stats_amino_acids,
        "get_protein_gene": get_protein_gene

    }
    if procedure == "check_mutations":
        results.append(check_mutations(seqs[0], seqs[1]))
    else:
        for seq in seqs:
            if is_protein(seq) is not True:
                raise ValueError("Invalid protein sequence")
            if procedure not in functions:
                raise ValueError("Wrong procedure name")
            else:
                results.append(functions[procedure](seq))
    if len(results) == 1:
        return results[0]
    else:
        return results
        