import itertools

def hamming_distance(s1, s2):
    """
    Calculate the Hamming distance between two strings of equal length.
    """
    if len(s1) != len(s2):
        raise ValueError("Input strings must have the same length")
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def levenshtein_distance(s1, s2):
    '''
    Levenshtein edit distance between 2 strings.
    Source: https://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Levenshtein_distance#Python
    '''
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1 # j+1 instead of j since previous_row and current_row are one character longer
            deletions = current_row[j] + 1       # than s2
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row
    return previous_row[-1]


def generate_hamming_strings(input_str, N, alphabet='ATCG'):
    """
    Generate all strings within a Hamming distance N of the input string.
    """
    result = set([input_str])
    if N <= 0:
        return result
    n = len(input_str)
    for i in range(n):
        for char in alphabet:
            modified_str = input_str[:i] + char + input_str[i + 1:]
            result |= generate_hamming_strings(modified_str, N - 1)
    return result


def generate_indel_strings(input_str, N, alphabet='ATCG'):
    """
    Generate all strings that can be obtained by introducing N indels (insertions or deletions)
    to the input string.
    """
    result = set([input_str])
    if N <= 0:
        return result

    # Insertion
    for i in range(len(input_str) + 1):
        for char in alphabet:
            new_str = input_str[:i] + char + input_str[i:]
            result |= generate_indel_strings(new_str, N - 1)

    # Deletion
    for i in range(len(input_str)):
        new_str = input_str[:i] + input_str[i + 1:]
        result |= generate_indel_strings(new_str, N - 1)

    return result


def generate_levenshtein_strings(seq, n_max, n_indel=None, n_subs=None, **kwargs):
    '''
    Generate set of strings within a specified edit distance from the input.

    Args
    - seq: str
        Input sequence
    - n_max: int
        Maximum Levenshtein edit distance
    - n_indel: int. default=None
        Maximum number of indels. If None, then n_index = n_max.
    - n_subs: int. default=None
        Maximum number of substitutions. If None, then n_subs = n_max.
    - **kwargs
        Additional arguments passed onto generate_hamming_strings() and generate_indel_strings()
        - alphabet: str

    Returns: set(str)
    '''
    assert n_max >= 0
    if n_indel is None:
        n_indel = n_max
    if n_subs is None:
        n_subs = n_max
    all_seqs = set([seq])
    for i in range(min(n_indel, n_max) + 1):
        hamming_seqs = set([seq])
        for j in range(min(n_subs, n_max - i) + 1):
            hamming_seqs |= generate_hamming_strings(seq, j, **kwargs)
        all_seqs |= hamming_seqs
        for hseq in hamming_seqs:
            all_seqs |= generate_indel_strings(hseq, i, **kwargs)
    return all_seqs


def generate_variant_map(seqs, dist_total, dist_hamming=None, dist_indel=None, verify_unique=True, variants_as_keys=True):
    '''
    Args
    - seqs: set(str)
        Sequences for which to generate variants
    - dist_total: int
        Maximum edit distance
    - dist_hamming: int. default=None
        Maximum Hamming distance
    - dist_indel: int. default=None
        Maximum number of indels
    - verify_unique: bool
        Assert that no 2 input sequences share any variants
    - variants_as_keys: bool. default=True
        Return a dict mapping from variants to their corresponding input sequence.
        If False, then map from input sequence to their variants.

    Returns: dict(str -> set(str))
        See variants_as_keys argument.
    '''
    result = {seq: generate_levenshtein_strings(seq, dist_total, n_indel=dist_indel, n_subs=dist_hamming) for seq in seqs}
    if verify_unique:
        current_seqs = set()
        for seq, variants in result.items():
            assert len(variants & current_seqs) == 0, \
                f'Non-unique variant encountered for sequence {seq}'
            current_seqs |= variants
    if variants_as_keys:
        return {variant: key for key, variants in result.items() for variant in variants}
    return result


def min_group_distance(seqs, distfun):
    '''
    Compute the minimum distance between any 2 sequences in a group.

    Args
    - seqs: iterable
        Sequences
    - distfun: callable
        Distance function, such as `hamming_distance` or `levenshtein_distance`.
        Must take 2 positional arguments (the 2 strings to be compared).

    Returns: int or np.nan
    - If distfun is `hamming_distance`, then np.nan is returned if not all
      the sequences in seqs have the same length.
    '''
    try:
        return min(distfun(a, b) for a, b in itertools.combinations(seqs, 2))
    except ValueError as e:
        if e.args[0] == 'Input strings must have the same length':
            return np.nan
        else:
            raise e