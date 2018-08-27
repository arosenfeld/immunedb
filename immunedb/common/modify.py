import csv
import dnautils
import sys

from immunedb.common.models import (Clone, NoResult, SampleMetadata, Sample,
                                    SampleStats, Sequence, SequenceCollapse)
from immunedb.identification.metadata import NA_VALUES
from immunedb.util.log import logger


def remove_duplicates(session, sample):
    logger.info('Removing duplicates from sample {}'.format(sample.id))
    all_seqs = session.query(Sequence).filter(
        Sequence.sample == sample
    ).order_by(
        Sequence.copy_number.desc()
    )

    buckets = {}
    for seq in all_seqs:
        key = (seq.v_gene, seq.j_gene, seq.cdr3_num_nts)
        buckets.setdefault(key, []).append(seq)

    for i, bucket in enumerate(buckets.values()):
        while len(bucket) > 0:
            larger = bucket.pop(0)
            for i in reversed(range(len(bucket))):
                smaller = bucket[i]
                if dnautils.equal(larger.sequence, smaller.sequence):
                    larger.copy_number += smaller.copy_number
                    session.delete(smaller)
                    del bucket[i]

    session.commit()


def update_metadata(session, args):
    SENTINEL = '__TEMP'  # Used to temporarily avoid duplicate name issues
    with open(args.new_metadata) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        new_meta = {l['name']: l for l in reader}

    # delete existing metadata
    sample_ids = {s.name: s.id for s in session.query(Sample).filter(
        Sample.name.in_(new_meta))}

    session.query(SampleMetadata).filter(
        SampleMetadata.sample_id.in_(sample_ids.values())
    ).delete(synchronize_session='fetch')

    ignore_fields = ['name', 'new_name', 'subject', 'file_name']
    for sample_name, row in new_meta.items():
        if sample_name not in sample_ids:
            logger.warning('No sample {} in database.  Ignoring.'.format(
                sample_name))
        sample_id = sample_ids[sample_name]
        logger.info('Updating metadata for {}'.format(row['name']))
        session.add_all([
            SampleMetadata(sample_id=sample_id, key=k, value=v)
            for k, v in row.items() if k not in ignore_fields and v not in
            NA_VALUES
        ])
        if row['new_name'] != row['name']:
            logger.info('  Updating sample name to {}'.format(row['new_name']))
            session.query(Sample).filter(Sample.name == row['name']).update({
                Sample.name: row['new_name'] + SENTINEL
            })

    logger.info('Verifying uniqueness')
    for sample in session.query(Sample).filter(
            Sample.name.like('%' + SENTINEL)):
        sample.name = sample.name[:-len(SENTINEL)]

    if session.query(Clone.id).filter(~Clone.tree.is_(None)).count() > 0:
        logger.warning('This database has at least one clonal lineage '
                       'constructed.  All lineages will need to be updated '
                       'to reflect the modified metadata.')
    session.commit()


def combine_samples(session, args):
    groups = {}

    for meta in session.query(SampleMetadata).filter(
            SampleMetadata.key == args.combine_field):
        groups.setdefault(meta.value, set()).add(meta.sample_id)
    all_subjects = set()
    for group_id, samples in groups.items():
        group_subs = session.query(Sample.subject_id).filter(
            Sample.id.in_(samples)
        ).group_by(Sample.subject_id)
        group_subs = [s.subject_id for s in group_subs]
        all_subjects.update(set(group_subs))
        if len(group_subs) > 1:
            logger.error('Cannot combine samples across subjects '
                         '(group "{}" has {} subjects)'.format(
                             group_id, len(group_subs)))
            sys.exit(1)

    all_samples = [s.id for s in session.query(Sample.id).filter(
        Sample.subject_id.in_(all_subjects))]

    logger.info('Resetting information for {} subjects'.format(
        len(all_subjects), len(all_samples)))
    logger.info('   Resetting collapsing')
    session.query(SequenceCollapse).filter(
        SequenceCollapse.sample_id.in_(all_samples)
    ).delete(synchronize_session=False)
    logger.info('   Resetting clones')
    session.query(Clone).filter(
        Clone.subject_id.in_(all_subjects)
    ).delete(synchronize_session=False)
    logger.info('   Resetting sample statistics')
    session.query(SampleStats).filter(
        SampleStats.sample_id.in_(all_samples)
    ).delete(synchronize_session=False)

    for group_id, samples in groups.items():
        final_sample_id = min(samples)
        logger.info('Combining {} samples into new sample "{}" (ID {})'.format(
            len(samples), group_id, final_sample_id))
        session.query(Sequence).filter(
            Sequence.sample_id.in_(samples)
        ).update({
            Sequence.sample_id: final_sample_id,
        }, synchronize_session=False)

        logger.info('Updating sample name and deleting empty samples')
        # collapse to one sample
        final_sample = session.query(Sample).get(final_sample_id)
        final_sample.name = group_id
        remove_duplicates(session, final_sample)

        logger.info('Moving noresults')
        session.query(NoResult).filter(
            NoResult.sample_id.in_(samples)
        ).update({
            'sample_id': final_sample_id
        }, synchronize_session=False)

        # delete the now-empty samples
        session.query(Sample).filter(
            Sample.id.in_(samples - set([final_sample_id]))
        ).delete(synchronize_session=False)

    session.commit()
    logger.info('Sequences successfully collapsed: please re-run '
                'immunedb_collapse and later pipeline steps.')
