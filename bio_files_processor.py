import os


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
                if name and seq: # every time a new name starts, we add full oneline sequence into the dict
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