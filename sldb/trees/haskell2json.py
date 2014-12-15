import json
import os
import shlex
import subprocess
from sldb.common.models import *

def traverse(root):
    data = root[0]

    muts = []
    for mut in data['nodeMutations']:
        pos, from_nt, to_nt = mut.split('-')
        muts.append({
            'pos': int(pos),
            'from': from_nt,
            'to': to_nt
        })

    tree = {
        'data': {
            'name': '',
            'mutations': muts,
            'sequence': data['nodeSequence'],
        },
        'children': []
    }
    children = root[1]
    for c in children:
        tree['children'].append(traverse(c))
    return tree


def run_haskell2json(session, args):
    '''Converts a haskell-generated JSON file into a more readable and usable
    tree for storage in the database'''
    if args.clone_ids is None:
        clones = map(lambda c: c.id, session.query(Clone.id).filter(Clone.tree
        is None))
    else:
        clones = args.clone_ids

    for clone in clones:
        print 'Processing clone {}'.format(clone)
        # Get the clone
        print 'Verifying clone'
        clone_inst = session.query(Clone).filter(Clone.id == clone).first()
        if clone is None:
            print '[ERROR] No clone with ID {} found'.format(clone)
            return

        # Path for FASTA file and haskell output tree
        fasta_path = '{}/clone_{}.fasta'.format(args.temp, clone)
        haskell_path = '{}/haskell_tree_{}.json'.format(args.temp, clone)

        print 'Writing FASTA file'
        # Create the FASTA file
        with open(fasta_path, 'w+') as fh:
            for seq in session.query(SequenceMapping.seq_id,
                                     SequenceMapping.sequence).filter(
                    SequenceMapping.clone_id == clone):
                fh.write('>{}\n{}\n'.format(seq.seq_id, seq.sequence))

        print 'Creating tree'
        # Run the tree program
        proc = subprocess.Popen(shlex.split('{} -i {} -o {}'.format(
            args.tree_prog, fasta_path, haskell_path)), stdout=subprocess.PIPE)
        stdout, stderr = proc.communicate()

        # We don't expect any errors
        if stderr is not None:
            print 'Got error from tree program:\n"{}"\nExiting'.format(stderr)
            os.remove(fasta_path)
            os.remove(haskell_path)
            return

        print 'Writing JSON'
        # Write the cleaned JSON tree to the database
        with open(haskell_path) as fh:
            clone_inst.tree = json.dumps(traverse(json.load(fh)))
        session.commit()
        print 'Cleaning up temporary files'
        os.remove(fasta_path)
        os.remove(haskell_path)
