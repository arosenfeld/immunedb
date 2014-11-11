import argparse
from collections import Counter

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine

from sldb.common.models import *
import sldb.util.lookups as lookups


def _consensus(strings):
    cons = []
    for chars in zip(*strings):
        cons.append(Counter(chars).most_common(1)[0][0])

    return ''.join(cons)


def _process_clones(session, per_commit, limit_subjects=None):
    q = session.query(Clone)
    if limit_subjects is not None:
        q = q.filter(Clone.subject_id.in_(limit_subjects))

    for i, clone in enumerate(q):
        seqs = session.query(Sequence).filter(
            Sequence.clone_id == clone.id).all()

        clone.cdr3_nt = _consensus(map(lambda e: e.junction_nt, seqs))
        clone.cdr3_aa = lookups.aas_from_nts(clone.cdr3_nt, '')

        if i > 0 and i % per_commit == 0:
            print 'Committed {}'.format(i)
            session.commit()

    session.commit()


def run_clones():
    parser = config.get_base_arg_parser('Generates consensus for clones')
    parser.add_argument('-c', type=int, default=1000, help='Number of'
                        ' clones to generate between commits')
    parser.add_argument('--subjects', nargs='+', type=int,
                        help='Limit generation to certain subjects')
    args = parser.parse_args()

    session = config.get_session(args)

    _process_clones(session, args.c, args.subjects)
