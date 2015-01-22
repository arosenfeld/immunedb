import re

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

def find_v_position(sequence, reverse=False):
    '''Tries to find the end of the V gene region'''
    if type(sequence) == str:
        sequence = Seq(sequence)
    # TODO: Possibly look for TATTACTG
    loc = sequence.rfind('TATTACTGT')
    if loc >= 0:
        return loc + 6
    # Try to find DxxxyzC
    found = _find_dc(sequence, reverse)
    if found is None:
        # If DxxyzC isn't found, try to find 'YYC', 'YCC', or 'YHC'
        found = _find_yxc(sequence, reverse)

    loc = sequence.rfind('TATTACTG')
    if loc >= 0:
        return loc + 6
    return found


def _find_dc(sequence, reverse):
    return _find_with_frameshifts(sequence, 'D(.{3}((YY)|(YC)|(YH)))C', reverse)


def _find_yxc(sequence, reverse):
    return _find_with_frameshifts(sequence, 'Y([YHC])C', reverse)


def _find_with_frameshifts(sequence, regex, reverse):
    r = range(2, -1, -1)
    if reverse:
        r = reversed(r)

    for shift in r:
        seq = sequence[shift:]
        seq = seq[:len(seq) - len(seq) % 3]
        aas = str(seq.translate())
        res = re.search(regex, aas)
        if res is not None:
            return (res.end() - 1) * 3 + shift
    return None
