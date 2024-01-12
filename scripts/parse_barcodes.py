import pandas as pd

def barcodes_to_df(f, regex, split='::', store_unmatched=100):
    '''
    Args
    - f: file object
        The barcode file, e.g., as the output of open() or gzip.open()
    - regex: re.Pattern
        Regular expression to search for in each line of f.
        Named groups specified with `(?P<name>...)` are extracted
        and used to generate the DataFrame.
    - split: str or None
        Split each line of the barcode file according to split,
        then the last split is used for regular expression searching.
    - store_unmatched: int. default=100
        The maximum number of unmatched lines to return. Useful for debugging.

    Returns
    - df: pd.DataFrame
        Column names are given by group names in regex.
    - n_unmatched: int
        Number of unmatched strings
    - unmatched: list of str
        Up to store_unmatched unmatched lines.
    '''
    barcodes = []
    n_unmatched = 0
    unmatched = []
    for line in f:
        target = line.strip().split(split)[-1]
        match = regex.search(target)
        if match is None:
            n_unmatched += 1
            if store_unmatched > 0:
                unmatched.append(line)
                store_unmatched -= 1
        else:
            barcodes.append(match.groupdict())
    return pd.DataFrame(barcodes), n_unmatched, unmatched

# def count_barcodes(df, rounds, col_umi='UMI'):
#     if type(col_umi) is not list:
        
#     col_umi = [col_umi] if type(col_umi) is not list else col_umi
#     df.groupby(rounds + col_umi).count()