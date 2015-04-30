import re

import distance

from sldb.common.models import Clone, CloneGroup, Sequence
import sldb.common.modification_log as mod_log
from sldb.identification.v_genes import VGene


def _update_clone(session, clone, v_name, new_v_seq):
    clone.v_gene = v_name
    group = session.query(CloneGroup).filter(
        CloneGroup.v_gene == v_name,
        CloneGroup.j_gene == clone.j_gene,
        CloneGroup.cdr3_aa == clone.group.cdr3_aa,
        CloneGroup.cdr3_num_nts == clone.cdr3_num_nts,
        CloneGroup.subject_id == clone.subject_id).first()
    if not group:
        germline = (new_v_seq[:VGene.CDR3_OFFSET:] +
                    clone.group.germline[VGene.CDR3_OFFSET:])
        group = CloneGroup(
            v_gene=v_name,
            j_gene=clone.j_gene,
            cdr3_aa=clone.group.cdr3_aa,
            cdr3_num_nts=clone.cdr3_num_nts,
            subject_id=clone.subject_id,
            germline=germline)
        session.add(group)

    # TODO: Update clone group
    clone.group = group


def _gap_to_germ(seq, germline):
    for i, c in enumerate(germline[:VGene.CDR3_OFFSET]):
        if c == '-':
            seq = seq[:i] + '-' + seq[i:]
    return seq


def _add_new_gaps(seq, gaps):
    added = 0
    for (pos, size) in gaps:
        seq = seq[:pos + added - 1] + ('-' * size) + seq[pos + added - 1:]
        added += size
    return seq


def _calculate_stats(seq):
    def _similarity(s1, s2):
        return len(s1) - sum([
            0 if gs == 'N' or gs == ss else 1 for gs, ss in zip(s1, s2)])

    seq.num_gaps = seq.sequence[0:VGene.CDR3_OFFSET].count('-')
    seq.pad_length = re.search('[ATCG]', seq.sequence).start()
    seq.v_match = _similarity(seq.sequence[seq.pad_length:VGene.CDR3_OFFSET],
                              seq.germline[seq.pad_length:VGene.CDR3_OFFSET])
    # TODO: This should count the streak into the CDR3
    seq.v_length = VGene.CDR3_OFFSET - seq.pad_length - seq.num_gaps
    seq.pre_cdr3_length = VGene.CDR3_OFFSET - seq.pad_length


def run_modify_clone(session, args):
    clone = session.query(Clone).filter(Clone.id == args.clone_id).first()
    if clone is None:
        print 'No clone with ID {}'.format(args.clone_id)
        return

    old_v = clone.v_gene
    old_v_seq = clone.group.germline
    if args.v_name:
        if clone.v_gene == args.v_name:
            print '[ERROR] V-gene specified is already assigned to clone'
            return
        existing = session.query(
            Sequence.sequence).filter(Sequence.v_gene == args.v_name).first()
        if existing:
            if args.v_seq and (existing.sequence[:VGene.CDR3_OFFSET] !=
                               args.v_seq[:VGene.CDR3_OFFSET].upper()):
                print '[ERROR] Existing V-gene specified with new sequence.'
                return
            else:
                print 'Using existing V gene and sequence {}'.format(
                    args.v_name)
                new_v_seq = existing.sequence[:VGene.CDR3_OFFSET]
        else:
            print 'Using new V gene and sequence {}'.format(args.v_name)
            new_v_seq = args.v_seq[:VGene.CDR3_OFFSET].upper()

        _update_clone(session, clone, args.v_name, new_v_seq)
        session.commit()

    if args.seq_gaps:
        gaps = []
        try:
            for gap in args.seq_gaps:
                pos, size = map(int, gap.split(',', 1))
                gaps.append((pos, size))
        except ValueError:
            print '[ERROR] Invalid gap string.'
            return
        seqs = session.query(Sequence).filter(Sequence.clone_id == clone.id)
        for i, seq in enumerate(seqs):
            if i > 0 and i % 1000 == 0:
                print 'Committed {}'.format(i)
                session.commit()
            new_seq = seq.sequence.replace('-', '')
            new_seq = _add_new_gaps(new_seq, sorted(gaps, key=lambda e: e[0]))
            new_seq = _gap_to_germ(new_seq, clone.group.germline)
            seq.sequence = new_seq
            seq.v_gene = clone.group.v_gene
            seq.germline = clone.group.germline
            _calculate_stats(seq)
        session.commit()
    session.add(mod_log.make_clone_mod(clone, old_v, old_v_seq, args.seq_gaps))
    session.commit()
