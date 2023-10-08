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
            if gc_content >= lower_bound and gc_content <= upper_bound:
                return True
            else:
                return False
    else:
        if gc_content < gc_bounds:
            return True
        else:
            return False


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
            if len(seq) >= lower_bound and len(seq) <= upper_bound:
                return True
            else:
                return False
    else:
        if len(seq) < length_bounds:
            return True
        else:
            return False


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
    
    if quality > quality_threshold:
        return True
    else:
        return False
