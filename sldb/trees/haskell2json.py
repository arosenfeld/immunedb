import json
import os
import shlex
import subprocess

from sqlalchemy.sql import func

from sldb.common.models import *


def traverse(session, root, start=False):
    data = root[0]

    muts = []
    for mut in data['nodeMutations']:
        pos, from_nt, to_nt = mut.split('-')
        muts.append({
            'pos': int(pos),
            'from': from_nt,
            'to': to_nt
        })

    seqs = filter(lambda d: d['remainingMutations'] == 0, data['sequences'])
    cn = 0
    for s in seqs:
        cn += int(s['printFastaHeader'].split('|')[1])

    if len(seqs) > 0:
        info = session.query(Sequence).join(Sample).filter(
            Sequence.seq_id.in_(map(
                lambda s: s['printFastaHeader'].split('|')[0], seqs)))
        info = info.all()
        tissues = ','.join(sorted(set(map(lambda s: s.sample.tissue, info))))
        subsets = ','.join(sorted(set(map(lambda s: s.sample.subset
                                          if s.sample.subset else '', info))))
    else:
        tissues = ''
        subsets = ''
    tree = {
        'data': {
            'mutations': muts,
            'copy_number': cn,
            'tissues': tissues,
            'subsets': subsets,
            'seq_ids': map(lambda s: s['printFastaHeader'].split('|')[0],
                           seqs),
            'sequence': data['nodeSequence'],
        },
        'children': []
    }
    children = root[1]
    for c in children:
        tree['children'].append(traverse(session, c))
    return tree


def run_haskell2json(session, args):
    '''Converts a haskell-generated JSON file into a more readable and usable
    tree for storage in the database'''
    if args.clone_ids is None:
        clones = session.query(Clone.id)
        if not args.force:
            clones = clones.filter(Clone.tree.is_(None))
        clones = map(lambda c: c.id, clones)
    else:
        clones = args.clone_ids

    for clone in clones:
        print 'Processing clone {}'.format(clone)
        # Get the clone
        clone_inst = session.query(Clone).filter(Clone.id == clone).first()
        if clone is None:
            print '[ERROR] No clone with ID {} found'.format(clone)
            return

        if clone_inst.tree is not None and not args.force:
            print ('Not regenerating tree for clone {}.  Use --force to '
                   'override.').format(clone)
           continue

        # Path for FASTA file and haskell output tree
        fasta_path = '{}/clone_{}.fasta'.format(args.temp, clone)
        haskell_path = '{}/haskell_tree_{}.json'.format(args.temp, clone)

        # Create the FASTA file
        with open(fasta_path, 'w+') as fh:
            fh.write('>{}\n{}\n'.format(
                'germline',
                session.query(CloneGroup.germline).filter(
                    CloneGroup.id == clone_inst.group_id).first().germline))
            for seq in session.query(
                    func.sum(Sequence.copy_number).label('copy_number'),
                    Sequence.seq_id,
                    Sequence.sequence,
                    Sequence.sequence_replaced).filter(
                        Sequence.clone_id == clone)\
                    .group_by(Sequence.sequence_replaced):
                fh.write('>{}|{}\n{}\n'.format(
                    seq.seq_id,
                    seq.copy_number,
                    seq.sequence_replaced))

        # Run the tree program
        proc = subprocess.Popen(shlex.split('{} -i {} -o {} -c -C 2'.format(
            args.tree_prog, fasta_path, haskell_path)), stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()

        # We don't expect any errors
        if stderr is not None:
            print 'Got error from tree program:\n"{}"\nExiting'.format(stderr)
            os.remove(fasta_path)
            os.remove(haskell_path)
            return

        # Write the cleaned JSON tree to the database
        with open(haskell_path) as fh:
            clone_inst.tree = json.dumps(traverse(session, json.load(fh)))
        session.commit()
        os.remove(fasta_path)
        os.remove(haskell_path)
