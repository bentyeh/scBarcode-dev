from collections import defaultdict
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import matplotlib.ticker
from helpers import fastq_parse

regex_loc_tag = re.compile(r"LX:Z:(([^:]+:\d+,\d+-\d+,?)+)")
regex_tag = re.compile(r"([^:]+):(\d+),(\d+)-(\d+),?")


def parse_locations(name):
    """
    Parse location tag from read name.

    Args
    - name: str
        Read name, including location tag.
        Example: '@readname LX:Z:tag_A:0,3-6,tag_B:0,6-10'

    Returns: list of 4-tuple (str, int, int, int), or None
    - If a properly formatted location tag is found, then returns a list where each element
      corresponds to a tag. Each tag is described by a 4-tuple: (tag name, file number,
      start position, end position), where positions are given as in Python indexing. Tags
      in the list are sorted by file number, start position, end position.
    - Otherwise returns None.
    """
    match = regex_loc_tag.search(name)
    if match is None:
        return None
    features = regex_tag.findall(match.group(1))
    features = sorted(
        [(name, int(file_n), int(start), int(end)) for (name, file_n, start, end) in features],
        key=lambda feature: (feature[1], feature[2], feature[3]),
    )
    return features


def features_to_coordinates(features):
    """
    Convert list of features to plotting coordinates

    Arg: list of 4-tuple (str, int, int, int)
    - Each tuple describes a feature: (name, file number, start position, end position)
    - Assumed to be sorted by file number, start position, and end position

    Returns: dict
    - Keys: file number
    - Values: list of 4-tuple (name, start, end, y)
    """
    name, file_n, start, end = features[0]
    current_y_ends = {file_n: [end]}
    feature_coords = defaultdict(list)
    feature_coords[file_n] = [(name, start, end, 0)]
    for name, file_n, start, end in features[1:]:
        if file_n in current_y_ends:
            use_current = False
            for y, y_end in enumerate(current_y_ends[file_n]):
                if start >= y_end:
                    use_current = True
                    current_y_ends[file_n][y] = end
                    break
            if not use_current:
                current_y_ends[file_n].append(end)
                y += 1
            feature_coords[file_n].append((name, start, end, y))
        else:
            current_y_ends[file_n] = [end]
            feature_coords[file_n].append((name, start, end, 0))
    return feature_coords


def parse_quals(qual_str, offset=33):
    """
    Args
    - qual_str: str
    - offset: int. default=33

    Returns: np.ndarray
    """
    return np.array([ord(c) - offset for c in qual_str])


def plot_features(feature_coords, ax):
    """
    Plot features from a read.

    Args
    - feature_coords: list of 4-tuple (str, int, int, int)
        Each tuple describes a feature: (name, start position, end position, y)
    - ax: matplotlib.axes.Axes

    Returns: list of matplotlib.patches.Rectangle, list of matplotlib.text.Text
    """
    rects = [
        Rectangle((start-0.5, y), (end - start), 1, facecolor=f"C{i}", edgecolor="black")
        for i, (name, start, end, y) in enumerate(feature_coords)
    ]
    for rect in rects:
        ax.add_patch(rect)
    texts = [
        ax.text((end - start) / 2 + start, y + 0.5, name, va="center", ha="center")
        for (name, start, end, y) in feature_coords
    ]
    ax.relim()
    ax.autoscale_view()
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.get_yaxis().set_visible(False)
    return rects, texts


def plot_seq_qscores(seq, quals, ax, offset=0):
    """
    Plot sequence and FASTQ scores.

    Args
    - seq: str
        Sequence
    - quals: np.ndarray
        Quality scores
    - ax: matplotlib.axes.Axes
    - offset: int. default=0
        Offset to use for x-axis ticks

    Returns: None
    """
    assert len(quals) == len(seq)
    xticks = range(offset, len(quals) + offset)
    ax.plot(xticks, quals)
    ax.set_xticks(xticks, labels=list(seq))
    ax.set_ylabel("Q score")


def plot_read(seq, quals, feature_coords, axs, reverse=False):
    """
    Plot sequence, FASTQ scores, and features from a read.

    Args
    - seq: str
        Sequence
    - quals: str
        Quality string
    - feature_coords: list of 4-tuple (str, int, int, int)
        Each tuple describes a feature: (name, start position, end position, y)
    - axs: list of matplotlib.axes.Axes
    - reverse: bool. default=False
        Reverse x-axis

    Returns: None
    """
    qscores = parse_quals(quals)
    plot_seq_qscores(seq, qscores, axs[0])
    plot_features(feature_coords, axs[1])
    xticks = list(range(len(seq)))
    axs[1].set_xticks(xticks, labels=xticks[::-1] if reverse else xticks)
    axs[1].xaxis.set_major_locator(matplotlib.ticker.AutoLocator())
    axs[0].margins(x=0.01)
    axs[1].set_xlim(*axs[0].get_xlim())


def plot_read_pair(name, seq1, quals1, seq2=None, quals2=None, reverse2=False, fig_kws=None):
    """
    Plot features from a read pair.

    Args
    - name: str
        Read name, including location tag.
        Example: '@readname LX:Z:tag_A:0,3-6,tag_B:0,6-10'
    - seq1: str
        Read 1 sequence
    - quals1: str
        Read 1 qualities
    - seq2: str. default=None
        Read 2 sequence
    - quals2: str
        Read 2 qualities
    - reverse2: bool. default=False
        Plot read 2 sequence, quality scores, and features in reverse orientation
        (but not reverse complement)
    - fig_kws: dict. default=None
        Keyword arguments to pass to plt.subplots().

    Returns: matplotlib.figure.Figure
    """
    if fig_kws is None:
        fig_kws = {}
    features = parse_locations(name)
    if features is None:
        # only plot sequence and FASTQ scores
        max_y = 0
    else:
        # plot sequence, FASTQ scores, features
        feature_coords = features_to_coordinates(features)
        max_y = max(
            [feature[3] + 1 for feature_list in feature_coords.values() for feature in feature_list]
        )
    fig_kws_default = dict(
        constrained_layout=True,
        height_ratios=[4, max_y / 2],
        width_ratios=[len(seq1), len(seq2)] if seq2 is not None else [1],
        sharey="row",
    )
    fig_kws_default.update(fig_kws)
    ncols = 1 if seq2 is None else 2
    fig, axs = plt.subplots(nrows=2, ncols=ncols, squeeze=False, **fig_kws_default)
    plot_read(seq1, quals1, feature_coords[0], [axs[0, 0], axs[1, 0]])
    axs[0, 0].set_title("read 1")
    if seq2:
        if reverse2:
            seq2 = seq2[::-1]
            quals2 = quals2[::-1]
            new_feature_coords = []
            n = len(seq2) - 1
            for feature_name, start, end, y in feature_coords[1]:
                start_rev = n - end
                end_rev = n - start
                new_feature_coords.append((feature_name, start_rev, end_rev, y))
            feature_coords[1] = new_feature_coords
        plot_read(seq2, quals2, feature_coords[1], [axs[0, 1], axs[1, 1]], reverse=reverse2)
        axs[0, 1].set_title("read 2")
    fig.suptitle(f"Read {name.split()[0]}")
    return fig
