from Bio.Seq import Seq

j_anchors = {
    'IGHJ1': 'TGGTCACCGTCTCCTCAG',
    'IGHJ2': 'TGGTCACTGTCTCCTCAG',
    'IGHJ3': 'TGGTCACCGTCTCTTCAG',
    'IGHJ4': 'TGGTCACCGTCTCCTCAG',
    'IGHJ5': 'TGGTCACCGTCTCCTCAG',
    'IGHJ6': 'CGGTCACCGTCTCCTCAG',
}

def all_j_anchors(min_length):
    max_size = max(map(len, j_anchors.values()))
    for trim in range(0, max_size, 3):
        for j, f in j_anchors.iteritems():
            new_f = f[:len(f) - trim]
            if len(new_f) >= min_length:
                yield new_f, f, j

def get_j_ties(j_name, j_seq):
    tied = set([j_name])
    for j, seq in sorted(j_anchors.iteritems()):
        if j == j_name:
            continue

        if seq[:len(j_seq)] == j_seq:
            tied.add(j)

    return tied
