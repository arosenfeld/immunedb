import dnautils
import re
from sldb.identification import AlignmentException
from sldb.identification.v_genes import get_common_seq
from sldb.common.models import CDR3_OFFSET

def _find_cdr3_start(seq):
    offset = 0
    for i, c in enumerate(seq):
        offset += 1 if c != '.' else 0
        if offset == CDR3_OFFSET:
            return i + 1

def align_v(sequence, v_germlines):
    max_align = None
    max_v = None
    for name, germr in v_germlines.iteritems():
        germ = germr.sequence
        try:
            g, s, score = dnautils.align(germ.replace('-', ''),
                                         sequence.replace('-', ''),
                                         -90, -30, -10, -10, 40)
            s = ''.join(reversed(s))
            g = ''.join(reversed(g))
        except Exception as e:
            continue
        if max_align is None or score >= max_align['score']:
            if max_align is None or score > max_align['score']:
                max_v = []
            max_v.append((name, g))
            if len(germ) > len(g):
                g = germ[:len(germ) - len(g)] + g
            max_align = {
                'original_germ': germ,
                'seq': s,
                'germ': g,
                'score': score
            }
    if max_align is None:
        raise AlignmentException('Could not locally align sequence.')

    germ_final = list(get_common_seq(
        map(lambda e: e[1], max_v), cutoff=False).replace('-', '.')
    )
    seq_final = list(max_align['seq'].replace('-', '.'))

    imgt_gaps = [
        i for i, c in enumerate(max_align['original_germ']) if c == '-'
    ]
    for gap in imgt_gaps:
        germ_final.insert(gap, '-')
        seq_final.insert(gap, '-')
    germ_final = ''.join(germ_final)
    seq_final = ''.join(seq_final)
    offset = re.search('[ATCGN]', germ_final)
    if offset is None:
        raise AlignmentException('Entire germline gapped.')
    offset = offset.start()
    germ_final = germ_final[offset:]
    seq_final = seq_final[offset:]

    end = _find_cdr3_start(germ_final)
    germ_final = germ_final[:end]#.replace('.', '-')
    seq_final = seq_final#.replace('.', '-')

    return map(lambda e: e[0], max_v), germ_final, seq_final
