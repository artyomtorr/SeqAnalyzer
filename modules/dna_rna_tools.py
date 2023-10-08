COMPLEMEMNT_DICT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'U': 'A',  
                 'a': 't', 't': 'a', 'g': 'c', 'c': 'g', 'u': 'a'}

TRANSCRIPTION_DICT = {'T': 'U', 't': 'u'}

ALPHABET_DNA = {'A', 'T', 'G', 'C'}

ALPHABET_RNA = {'A', 'G', 'C', 'U'}


def reverse(seq:str) -> str:
    '''
    Returns reversed nucleotide sequence.
    '''
    return seq[::-1]


def complement(seq:str) ->str:
    '''
    Returns complementary nucleotide sequence.
    '''
    result = []
    for nucleotide in seq:
        if nucleotide in COMPLEMEMNT_DICT:
            result += COMPLEMEMNT_DICT[nucleotide]
        else: 
            result += nucleotide
    return "".join(result)


def reverse_complement(seq:str) ->str:
    '''
    Returns reversed complementary nucleotide sequence.
    '''
    return complement(reverse(seq))


def transcribe(seq:str)-> str:
    '''
    Returns transcribed nucleotide sequence.
    '''
    result = []
    for nucleotide in seq:
        if nucleotide in TRANSCRIPTION_DICT:
            result += TRANSCRIPTION_DICT[nucleotide]
        else: 
            result += nucleotide
    return "".join(result)


def is_rna(seq:str) -> bool:
    """
    Check if the RNA sequence consist only RNA nucleotides.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= ALPHABET_RNA


def is_dna(seq:str) -> bool:
    """
    Check if the DNA sequence consist only DNA nucleotides.
    """
    unique_chars = set(seq.upper())
    return unique_chars <= ALPHABET_DNA
