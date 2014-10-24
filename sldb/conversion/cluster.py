import argparse
from collections import Counter

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine, func
from sqlalchemy.engine.reflection import Inspector

from sldb.common.models import *


def _similar(s1, s2, min_similarity):
    assert len(s1) == len(s2)
    return sum(ch1 == ch2 for ch1, ch2 in zip(s1, s2)) / len(s1) >= \
        min_similarity


def _get_groups_from_db(session):
    for g in session.query(Sequence.v_call, Sequence.j_call,
                         func.count(Sequence.junction_aa)
                         .label('cdr3_len'))\
        .group_by(Sequence.v_call, Sequence.j_call,
                  func.count(Sequence.junction_aa)):
        yield (g.v_call, g.j_call, g.cdr3_len)


def _get_groups_from_file(fn):
    with open(fn) as fh:
        for l in fh:
            fields = l.split('\t')
            yield (fields[0], fields[1], int(fields[2]))


def _intragroup_similar(seqs, pseq, min_similarity):
    for seq in seqs:
        if not _similar(pseq, seq, min_similarity):
            return False

    return True


def _cluster_group(session, group, min_similarity):
    seq_clusters = []
    for seq in session.query(Sequence.junction_aa)\
        .filter(Sequence.v_call == group[0])\
        .filter(Sequence.j_call == group[1])\
        .filter(func.length(Sequence.junction_nt) == group[2]).yield_per(100000):

        cdr3 = seq.junction_aa

        added = False
        for seqs in seq_clusters:
            if _similar(next(iter(seqs)), cdr3, min_similarity):
                if _intragroup_similar(seqs, cdr3, min_similarity):
                    seqs.add(cdr3)
                    added = True
                    break

        if not added:
            seq_clusters.append(set([cdr3]))

    return seq_clusters


def _generate_consensus_cdr3(cluster):
    cons = []
    for chars in zip(*cluster):
        cons.append(Counter(chars).most_common(1)[0][0])

    return ''.join(cons)


def run_cluster():
    parser = argparse.ArgumentParser(description='Runs clonal identification')
    parser.add_argument('host', help='mySQL host')
    parser.add_argument('db', help='mySQL database')
    parser.add_argument('user', help='mySQL user')
    parser.add_argument('pw', help='mySQL password')
    parser.add_argument('res_dir', help='Directory for cluster output')
    parser.add_argument('-p', type=int, default=85, help='Percent overlap to '
                        'consider two CDR3s as part of the same clone')
    parser.add_argument('-i', default=None, help='If specified will use the '
                        'parameter as groupings instead of querying database')
    args = parser.parse_args()

    engine = create_engine(('mysql://{}:{}@{}/'
                            '{}?charset=utf8&use_unicode=0').format(
                                args.user, args.pw, args.host, args.db))

    Base.metadata.create_all(engine)
    Base.metadata.bind = engine
    session = sessionmaker(bind=engine)()

    if args.i == None:
        groups = _get_groups_from_db(session)
    else:
        groups = _get_groups_from_file(args.i)

    for group in groups:
        v = group[0]
        j = group[1]
        cdr3_len = group[2]
        print 'Clustering {} {} {}'.format(v, j, cdr3_len)
        clusters = _cluster_group(session, group, args.p / 100.0)

        path = '{}/{}__{}__{}.cluster'.format(
            args.res_dir,
            v.replace('/','_'),
            j.replace('/','_'),
            cdr3_len)

        with open(path, 'w+') as fh:
            for c in clusters:
                fh.write('{}\t{}\n'.format(
                    _generate_consensus_cdr3(c), ','.join(c)))
