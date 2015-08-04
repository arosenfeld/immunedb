import dnautils
import re

from sqlalchemy import func

from sldb.identification import AlignmentException
from sldb.identification.identify import SequenceRecord
from sldb.identification.j_genes import JGermlines
from sldb.identification.v_genes import get_common_seq, VGermlines
from sldb.identification.vdj_sequence import VDJSequence
from sldb.common.models import CDR3_OFFSET, DuplicateSequence, Sequence


def _find_cdr3_start(seq):
    offset = 0
    for i, c in enumerate(seq):
        offset += 1 if c != '.' else 0
        if offset == CDR3_OFFSET:
            return i + 1


def align_v(sequence, v_germlines, insert_penalty=-30, delete_penalty=-30,
            extend_penalty=-10, mismatch_penalty=-10, match_score=40):
    max_align = None
    max_v = None
    for name, germr in v_germlines.iteritems():
        germ = germr.sequence
        try:
            g, s, germ_omit, seq_omit, score = dnautils.align(
                germ.replace('-', ''), sequence.replace('-', ''),
                insert_penalty, delete_penalty, extend_penalty,
                mismatch_penalty, match_score)
            s = ''.join(reversed(s))
            g = ''.join(reversed(g))
        except Exception as e:
            continue
        if max_align is None or score >= max_align['score']:
            if max_align is None or score > max_align['score']:
                max_v = []
            max_v.append((name, g))
            if len(germ) > len(g):
                g = germ[:len(germ) - len(g) - germ_omit] + g
            if len(sequence) > len(s):
                s += sequence[-(len(sequence) - seq_omit - len(
                    s.replace('-', '').replace('.', ''))):]
            max_align = {
                'original_germ': germ,
                'seq': s,
                'germ': g,
                'score': score,
                'germ_omit': germ_omit,
                'seq_omit': seq_omit,
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
    germ_final = germ_final[:end]
    seq_final = seq_final

    return map(lambda e: e[0], max_v), germ_final, seq_final, end


def run_fix_indels(session, args):
    v_germlines = VGermlines(args.v_germlines)
    j_germlines = JGermlines(args.j_germlines, args.upstream_of_cdr3,
                             args.anchor_len, args.min_anchor_len)

    indels = session.query(Sequence).filter(
        Sequence.probable_indel_or_misalign == 1)
    total = indels.count()
    fixed = 0

    for i, seq in enumerate(indels):
        quality = seq.quality
        if quality is not None:
            quality = seq.quality.replace(' ', '')
            quality = map(lambda q: ord(q) - 33, quality)

        v = VDJSequence(seq.seq_id, seq.sequence.replace('-', ''), v_germlines,
                        j_germlines, quality, skip_align=True)
        try:
            avg_mut, avg_len = session.query(
                1 - func.avg(Sequence.v_match / Sequence.v_length),
                func.avg(Sequence.v_length)
            ).filter(Sequence.sample == seq.sample).first()
            v.locally_align(align_v(seq.sequence.replace('-', ''),
                                    v_germlines), avg_mut, avg_len)
            if not v.has_possible_indel:
                fixed += 1
                r = SequenceRecord(seq.sequence.replace('-', ''), seq.sample)
                r.seq_ids = map(lambda e: e.seq_id, session.query(
                    DuplicateSequence.seq_id
                ).filter(
                    DuplicateSequence.duplicate_seq_id == seq.seq_id
                ).all())
                r.seq_ids.append(v.id)
                r.vdj = v
                session.delete(seq)
                r.add_as_sequence(session, seq.sample, seq.paired)
        except AlignmentException as e:
            pass
        if i > 0 and i % 10 == 0:
            print 'Aligned {}/{} (fixed {}, {}%)'.format(
                i, total, fixed, round(100 * fixed / i))
    session.commit()
