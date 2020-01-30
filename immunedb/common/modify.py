import csv
import dnautils
import sys

from immunedb.common.models import (Clone, NoResult, SampleMetadata, Sample,
                                    Sequence, Subject)
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
        key = (seq.v_gene, seq.j_gene, seq.cdr3_num_nts, seq._insertions,
               seq._deletions)
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
    IGNORE_FIELDS = ['id', 'name', 'subject']
    with open(args.new_metadata) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        new_meta = {int(l['id']): l for l in reader}

    session.query(SampleMetadata).filter(
        SampleMetadata.sample_id.in_(new_meta.keys())
    ).delete(synchronize_session='fetch')

    old_subjects = set()
    for sample_id, row in new_meta.items():
        logger.info('Updating metadata for #{id}: {name}'.format(
            **row))
        sample = session.query(Sample).get(sample_id)
        # Update subject
        if sample.subject.identifier != row['subject']:
            logger.info('Subject for sample "{}" changed from {} -> {}'.format(
                sample.name, sample.subject.identifier, row['subject']
            ))
            old_subjects.add(sample.subject)
            new_subject = session.query(Subject).filter(
                Subject.study == sample.study,
                Subject.identifier == row['subject']
            ).first()
            if not new_subject:
                new_subject = Subject(
                    study=sample.study,
                    identifier=row['subject']
                )
                session.add(new_subject)
                session.flush()
                logger.info('\tNew subject found')
            else:
                old_subjects.add(new_subject)

            sample.subject = new_subject
            assert new_subject.id is not None
            session.query(Sequence).filter(Sequence.sample == sample).update({
                'subject_id': new_subject.id
            })

        for subject in session.query(Subject):
            if not subject.samples:
                logger.info('Deleting orphan subject "{}"'.format(
                    subject.identifier))
                session.delete(subject)
            elif subject in old_subjects:
                logger.info('Resetting subject "{}"'.format(
                    subject.identifier))
                subject.reset()

        # Update metadata
        session.add_all([
            SampleMetadata(sample=sample, key=k, value=v)
            for k, v in row.items() if k not in IGNORE_FIELDS and v not in
            NA_VALUES
        ])
        # Update name
        sample.name = row['name'] + SENTINEL

    session.commit()

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

    subjects = set()
    for meta in session.query(SampleMetadata).filter(
            SampleMetadata.key == args.combine_field):
        groups.setdefault(meta.value, set()).add(meta.sample)
        subjects.add(meta.sample.subject)

    for group_id, samples in groups.items():
        group_subs = set(s.subject for s in samples)
        if len(group_subs) > 1:
            logger.error('Cannot combine samples across subjects '
                         '(group "{}" has {} subjects)'.format(
                             group_id, len(group_subs)))
            sys.exit(1)

    for subject in subjects:
        subject.reset()

    for group_id, samples in groups.items():
        all_sample_ids = set(s.id for s in samples)
        final_sample_id = min(all_sample_ids)
        logger.info('Combining {} samples into new sample "{}" (ID {})'.format(
            len(samples), group_id, final_sample_id))
        session.query(Sequence).filter(
            Sequence.sample_id.in_(all_sample_ids)
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
            NoResult.sample_id.in_(all_sample_ids)
        ).update({
            'sample_id': final_sample_id
        }, synchronize_session=False)

        # delete the now-empty samples
        session.query(Sample).filter(
            Sample.id.in_(all_sample_ids - set([final_sample_id]))
        ).delete(synchronize_session=False)

    session.commit()
    logger.info('Sequences successfully collapsed: please re-run '
                'immunedb_collapse and later pipeline steps.')


def delete_samples(session, args):
    for sample_id in args.sample_ids:
        sample = session.query(Sample).get(sample_id)
        if not sample:
            logger.warning('Sample #{} does not exist'.format(sample_id))
            continue

        logger.info('Deleting sample #{}'.format(sample_id))
        sample.subject.reset()
        session.query(Sequence).filter(Sequence.sample == sample).delete()
        session.delete(sample)
    session.commit()

    for subject in session.query(Subject):
        if not subject.samples:
            logger.info('Deleting orphan subject "{}"'.format(
                subject.identifier))
            session.delete(subject)
    session.commit()
