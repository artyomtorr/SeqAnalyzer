import modules.fastq_filter as ff
import modules.dna_rna_tools as drt
import modules.protein_tools as pt


def run_fastq_filter(input_path: str, 
                 output_filename: str = None, 
                 gc_bounds: Tuple[int, int] = (0, 100), 
                 length_bounds: Tuple[int, int] = (0, 2**32), 
                 quality_threshold: int = 0) -> None:
    '''
    Filter FASTQ-sequences based on entered requirements.
    
    Arguments:
        - input_path (str): path to the file with FASTQ-sequences
        - output_filename (str): name of the output file with 
        filtered FASTQ-sequences
        - gc_bounds (tuple or int, default = (0, 100)): GC-content
        interval (percentage) for filtering. Tuple if contains 
        lower and upper bounds, int if only contains an upper bound.
        - length_bounds (tuple or int, default = (0, 2**32)): length 
        interval for filtering. Tuple if contains lower and upper 
        bounds, int if only contains an upper bound.
        - quality_threshold (int, default = 0): threshold value of average 
    read quality for filtering.

    Note: the output file is saved to the /fastq_filtrator_results 
    directory. The default output file name is the name of the input file.
    '''
    seqs = ff.read_fastq_file(input_path)
    filtered_seqs = {}
    for name in seqs:
        quality_value = ff.check_quality(seqs[name][2], quality_threshold) 
        length_value = ff.check_length(seqs[name][0], length_bounds)
        gc_content_value = ff.check_gc_content(seqs[name][0], gc_bounds)
        if  quality_value and length_value and gc_content_value:
            filtered_seqs[name] = seqs[name]
            
    if not output_filename:
        output_filename_with_extension = os.path.basename(input_path) 
        output_filename = os.path.splitext(output_filename_with_extension)[0] 
    ff.write_fastq_file(filtered_seqs, output_filename)


def run_dna_rna_tools(*args: str):
    """
    Function containing methods for DNA or RNA sequences analysis.

    Takes arbitrary number of arguments with DNA or RNA sequencies
    and the name of the procedure to be performed (always the last
    argument). 
    Returns the result of the procedure as string if one sequence
    is submitted or list if several.

    Supported procedures:
    - 'transcribe' - transcribes DNA sequence to RNA sequence.
    - 'reverse' - reverses nucleotide sequence.
    - 'complement' - returns complementary nucleotide sequence.
    - 'reverse_complement' - returns reversed complementary nucleotide sequence.
    """
    *seqs, procedure = args
    results = []
    procedures = {'transcribe': drt.transcribe, 
                 'reverse': drt.reverse,
                 'reverse_complement': drt.reverse_complement,
                 'complement': drt.complement}
    for seq in seqs:
        if not is_rna(seq) and is_dna(seq):
            raise ValueError("Invalid nucleotide sequence")
        results.append(procedures[procedure](seq))
    if len(results) == 1:
        return results[0]
    else: 
        return results


def run_protein_tools(*args: str):
    """
    Function containing methods for protein analysis.

    Takes arbitrary number of arguments with protein sequencies
    and the name of the procedure to be performed (always the last
    argument). Returns the result of the procedure as string
    or dictionary if one sequence is submitted or list if several.

    Note: if procedure 'check_mutations' is used then input must
    contain only three arguments: RNA sequence, protein sequence
    and the name of procedure itself.

    Supported procedures:
    - 'compute_molecular_weight' — computes molecular weight in g/mol
    - 'compute_length' — computes the number of amino acids
    - 'compute_hydrophobicity' — computes the percentage of gydrophobic aminoacids 
    - 'check_mutations' — checks missense mutations after translation
    - 'get_protein_gene' - returns possible variants of coding DNA sequences
    - 'count_amino_acids' - calculates the number of each aminoacid
    """
    *seqs, procedure = args
    results = []
    procedures = {
        "compute_molecular_weight": pt.compute_molecular_weight,
        "compute_length": pt.compute_length,
        "compute_hydrophobicity": pt.compute_hydrophobicity,
        "stats_amino_acids": pt.stats_amino_acids,
        "get_protein_gene": pt.get_protein_gene

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
                results.append(procedures[procedure](seq))
    if len(results) == 1:
        return results[0]
    else:
        return results        
