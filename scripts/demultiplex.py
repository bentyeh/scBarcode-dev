import sys
sys.path.append('../')
import Bio
import re
import string_distances

regex_Ns = re.compile('N+', flags=re.IGNORECASE)

def get_aligned_target_coordinates(alignment):
    return (alignment.coordinates[0, 0], alignment.coordinates[0, -1])

def index_alignments(adapters, aligner=None, regex_index=None):
    '''
    Alignments denoting where an index is in an adapter.

    Args
    - adapters: dict(str -> str)
        Map from name of adapter to adapter sequence
    - aligner: Bio.Align.PairwiseAligner. default=None
        Aligner

    Returns: dict(str -> Bio.Align.Alignment)
      Map from adapter name to alignment of index pattern (query) to adapter (target)
    '''
    if aligner is None:
        aligner = Bio.Align.PairwiseAligner(mismatch_score=-1, internal_gap_score=-1)
    if regex_index is None:
        regex_index = regex_Ns
    index_alignments = {}
    for name, seq in adapters.items():
        Ns = regex_index.search(seq).group()
        index_alignments[name] = aligner.align(seq, Ns)[0]
    return index_alignments


def extract_index(adapter_alignments, index_alignments, indices_hash=None, sort=True):
    '''
    Args
    - adapter_alignments: sequence((adapter name, Bio.Align.Alignment))
        Sequence of tuples of adapter name and adapter alignment (e.g., to a read)
    - index_alignments: dict(str -> Bio.Align.Alignment)
        Map from adapter name to alignment of index pattern (query) to adapter (target).
    - indices_hash: dict(str -> str). default=None
        Map from index sequence to assigned index name.
        If not provided, the index sequence is used as the index name.
    - sort: bool. default=True
        Sort returned alignments by target coordinate.

    Returns: list(tuple(int, int, str, str, str))
      List of adapter matches and associated indices
      - start coordinate
      - end coordinate
      - adapter name
      - index label
      - index sequence
    '''
    results = []
    for name, alignment in adapter_alignments:
        index_seq = alignment.map(index_alignments[name])[0]
        if indices_hash is not None:
            index_label = indices_hash[index_seq]
        else:
            index_label = index_seq
        results.append((alignment.coordinates[0, 0], alignment.coordinates[0, -1], name, index_label, index_seq))
    if sort:
        results.sort()
    return results


def find_adapters(
    read,
    adapters,
    thresholds,
    aligner=None,
    collapse_identical_coordinates=False,
    find_all_alignments_above_threshold=False):
    '''
    Args
    - read: str
    - adapters: dict(str -> str)
        Map from adapter name to adapter sequence, with a string of Ns denoting adapter indices
    - thresholds: dict(str -> numeric)
        Map from adapter name to adapter alignment score threshold (>=)
    - aligner: Bio.Align.PairwiseAligner
    - collapse_identical_coordinates: bool. default=False

    Returns: sequence((adapter name, Bio.Align.Alignment))
        Sequence of tuples of adapter name and adapter alignment (e.g., to a read)

    Caveat: currently only returns alignments with the same best score
    (a limitation of Bio.Align.PairwiseAligner), so if there are 2 different
    locations in a read where an adapter can match, this currently will only find
    one of them.
    - Potential solutions
      - Modify Bio.Align.PairwiseAligner to alignments with scores above a threshold
      - Perform global alignment to sliding windows across the read
      - After the first alignment, perform alignment to sequences on either side of the first alignment
    '''
    if aligner is None:
        aligner = Bio.Align.PairwiseAligner(mismatch_score=-1, internal_gap_score=-1, wildcard='N')
    adapter_alignments = []
    for name, adapter_seq in adapters:
        alignments = aligner.align(read, adapter_seq)
        if alignments.score >= thresholds[name]:
            if collapse_identical_coordinates:
                alignments.sort(key=get_aligned_target_coordinates)
                adapter_alignments.append(alignments[0])
                current_coords = get_aligned_target_coordinates(alignments[0])
                for alignment in alignments[1:]:
                    coords = get_aligned_target_coordinates(alignment)
                    if coords == current_coords:
                        continue
                    else:
                        current_coords = coords
                        adapter_alignments.append((name, alignment))
            else:
                for alignment in alignments:
                    adapter_alignments.append((name, alignment))
            # if find_all_alignments_above_threshold:
            #     aligned_coordinates = np.vstack(list(set(map(get_aligned_target_coordinates, alignments))))
            #     aligned_coordinates_outer = aligned_coordinates[:, 0].min()
            #     unaligned_segments = set(map(get_aligned_target_coordinates, alignments))
            #     find_adapters(
            #         read,
            #         adapters,
            #         thresholds,
            #         aligner=aligner,
            #         collapse_identical_coordinates=collapse_identical_coordinates,
            #         find_all_alignments_above_threshold=True)
    return adapter_alignments


def demultiplex(record, mod_names=False, loc_names=False, file=True, mapping=None):
    alignments = extract_index()
    if loc_names:
        pass
    if mod_names:
        pass
    if mapping:
        pass
    if file is True:
        filename = ''
        for (start, end, adapter_name, index_label, index_seq) in alignments:
            filename += f'({adapter_name}-{index_label})'
        with open(filename, 'a') as f:
            f.write()
