import datetime

from sqlalchemy import (Column, Boolean, Integer, String, Text, Date, DateTime,
                        ForeignKey, UniqueConstraint, Index, func)
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.orm.interfaces import MapperExtension
from sqlalchemy.dialects.mysql import TEXT, MEDIUMTEXT

import sldb.util.funcs as funcs
from sldb.common.settings import DATABASE_SETTINGS

if DATABASE_SETTINGS['master_metadata'] is None:
    print 'WARNING: Most master metadata specified'
if DATABASE_SETTINGS['data_metadata'] is None:
    print 'WARNING: Most data metadata specified'

BaseMaster = declarative_base(
    metadata=DATABASE_SETTINGS['master_metadata'])
BaseData = declarative_base(
    metadata=DATABASE_SETTINGS['data_metadata'])


class Study(BaseMaster):
    """A high-level study such as one studying Lupus.

    :param int id: An auto-assigned unique identifier for the study
    :param str name: A unique name for the study
    :param str info: Optional information about the study

    """
    __tablename__ = 'studies'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    id = Column(Integer, primary_key=True)
    # The name of the study
    name = Column(String(length=128), unique=True)
    # Some arbitrary information if necessary
    info = Column(String(length=1024))


class Subject(BaseMaster):
    """A subject which was sampled for a study.

    :param int id: An auto-assigned unique identifier for the subject
    :param str identifier: An identifier for the subject as defined by \
        the experimenter
    :param int study_id: The ID of the study under which the subject was \
        sampled
    :param Relationship study: Reference to the associated :py:class:`Study` \
        instance

    """
    __tablename__ = 'subjects'
    __table_args__ = (UniqueConstraint('study_id', 'identifier'),
                      {'mysql_engine': 'TokuDB'})

    id = Column(Integer, primary_key=True)

    identifier = Column(String(64))
    study_id = Column(Integer, ForeignKey(Study.id))
    study = relationship(Study, backref=backref('subjects',
                         order_by=identifier))


class Sample(BaseMaster):
    """A sample taken from a single subject, tissue, and subset.

    :param int id: An auto-assigned unique identifier for the sample
    :param str name: A unique name for the sample as defined by the \
        experimenter
    :param str info: Optional information about the sample
    :param date date: The date the sample was taken

    :param int study_id: The ID of the study under which the subject was \
        sampled
    :param Relationship study: Reference to the associated :py:class:`Study` \
        instance

    :param int subject_id: The ID of the subject from which the sample was \
        taken
    :param Relationship subject: Reference to the associated \
        :py:class:`Subject` instance


    :param str tissue: The tissue of the sample
    :param str subset: The tissue subset of the sample
    :param str disease: The known disease(s) present in the sample
    :param str lab: The lab which acquired the sample
    :param str experimenter: The experimenters name who took the sample

    """
    __tablename__ = 'samples'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    id = Column(Integer, primary_key=True)
    name = Column(String(128), unique=True)
    info = Column(String(1024))

    date = Column(Date)

    study_id = Column(Integer, ForeignKey(Study.id))
    study = relationship(Study, backref=backref('samples', order_by=(date,
                                                name)))

    subject_id = Column(Integer, ForeignKey(Subject.id), index=True)
    subject = relationship(Subject, backref=backref('samples',
                           order_by=(id)))
    subset = Column(String(128))
    tissue = Column(String(16))
    disease = Column(String(32))
    lab = Column(String(128))
    experimenter = Column(String(128))


class SampleStats(BaseData):
    """Aggregate statistics for a sample.  This exists to reduce the time
    queries take for a sample.

    :param int sample_id: The ID of the sample for which the statistics were \
        generated
    :param Relationship sample: Reference to the associated \
        :py:class:`Sample` instance

    :param str filter_type: The type of filter for the statistics
        (e.g. functional)
    :param bool outliers: If outliers were included in the statistics
    :param bool full_reads: If only full reads were included in the statistics
    :param str v_match_dist: Distribution of V gene match count
    :param str v_length_dist: Distribution of V gene total length
    :param str j_match_dist: Distribution of J gene match count
    :param str j_length_dist: Distribution of J gene total length
    :param str v_gene_dist: Distribution of V-gene assignments
    :param str j_gene_dist: Distribution of J-gene assignments
    :param str copy_number_dist: Distribution of copy number
    :param str cdr3_length_dist: Distribution of CDR3 lengths

    :param int sequence_cnt: The total number of sequences
    :param int in_frame_cnt: The number of in-frame sequences
    :param int stop_cnt: The number of sequences containing a stop codon
    :param int functional_cnt: The number of functional sequences
    :param int no_result_cnt: The number of invalid sequences

    """
    __tablename__ = 'sample_stats'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    sample_id = Column(Integer, ForeignKey(Sample.id),
                       primary_key=True)
    sample = relationship(Sample, backref=backref('sample_stats',
                          order_by=sample_id))

    filter_type = Column(String(length=255), primary_key=True)
    outliers = Column(Boolean, primary_key=True)
    full_reads = Column(Boolean, primary_key=True)

    v_identity_dist = Column(MEDIUMTEXT)

    v_match_dist = Column(MEDIUMTEXT)
    v_length_dist = Column(MEDIUMTEXT)

    j_match_dist = Column(MEDIUMTEXT)
    j_length_dist = Column(MEDIUMTEXT)

    v_gene_dist = Column(MEDIUMTEXT)
    j_gene_dist = Column(MEDIUMTEXT)

    copy_number_dist = Column(MEDIUMTEXT)
    cdr3_length_dist = Column(MEDIUMTEXT)

    sequence_cnt = Column(Integer)
    in_frame_cnt = Column(Integer)
    stop_cnt = Column(Integer)
    functional_cnt = Column(Integer)
    no_result_cnt = Column(Integer)


class CloneGroup(BaseMaster):
    """A group of clones which share the same V, J, and CDR3 amino-acids.  This
    is used to correlate identical or similar clones across database versions
    when their ID may change.

    :param int id: An auto-assigned unique identifier for the group
    :param str v_gene: The V-gene assigned to the sequence
    :param str j_gene: The J-gene assigned to the sequence
    :param cdr3_aa: The amino-acid sequence of the group's CDR3
    :param cdr3_num_nts: The number of nucleotides in the group's CDR3
    :param int subject_id: The ID of the subject from which the sample was \
        taken
    :param Relationship subject: Reference to the associated \
        :py:class:`Subject` instance
    :param str germline: The germline sequence for this sequence

    """
    __tablename__ = 'clone_groups'
    __table_args__ = (Index('grp_aas', 'v_gene', 'j_gene', 'subject_id',
                            'cdr3_aa'),
                      Index('grp_len', 'v_gene', 'j_gene', 'subject_id',
                            'cdr3_num_nts'),
                      UniqueConstraint('v_gene', 'j_gene', 'subject_id',
                                       'cdr3_num_nts', 'cdr3_aa'),
                      {'mysql_engine': 'TokuDB'})

    id = Column(Integer, primary_key=True)

    v_gene = Column(String(length=512))
    j_gene = Column(String(length=128))
    cdr3_aa = Column(String(length=128))
    cdr3_num_nts = Column(Integer)

    subject_id = Column(Integer, ForeignKey(Subject.id), index=True)
    subject = relationship(Subject, backref=backref('clones',
                           order_by=(v_gene, j_gene, cdr3_num_nts, cdr3_aa)))

    germline = Column(String(length=1024))


class Clone(BaseData):
    """A group of sequences likely originating from the same germline

    :param int id: An auto-assigned unique identifier for the clone
    :param str cdr3_nt: The consensus nucleotides for the clone
    :param int group_id: The group ID for the clone to correlate clones \
        across database versions

    :param Relationship group: Reference to the associated \
        :py:class:`CloneGroup` instance
    :param str tree: The textual representation of the clone's lineage tree

    """
    __tablename__ = 'clones'
    __table_args__ = (Index('bucket', 'v_gene', 'j_gene',
                            'cdr3_num_nts', 'subject_id'),
                      {'mysql_engine': 'TokuDB'})
    id = Column(Integer, primary_key=True)
    cdr3_nt = Column(String(length=512))

    # These are necessary during creation of clones, but will
    # always be redundant with the associated grouping after completion
    v_gene = Column(String(length=512), index=True)
    j_gene = Column(String(length=128), index=True)
    cdr3_num_nts = Column(Integer, index=True)
    subject_id = Column(Integer, ForeignKey(Subject.id), index=True)
    #

    group_id = Column(Integer, ForeignKey(CloneGroup.id),
                      index=True)
    group = relationship(CloneGroup, backref=backref('clones',
                         order_by=(id)))
    tree = Column(MEDIUMTEXT)


class CloneStats(BaseData):
    __tablename__ = 'clone_stats'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    clone_id = Column(Integer, ForeignKey(Clone.id), primary_key=True)
    clone = relationship(Clone)

    sample_id = Column(Integer, ForeignKey(Sample.id), primary_key=True)
    sample = relationship(Sample, backref=backref('clone_stats'))

    unique_cnt = Column(Integer)
    total_cnt = Column(Integer)

    mutations = Column(MEDIUMTEXT)
    selection_pressure = Column(MEDIUMTEXT)


class SequenceExtension(MapperExtension):
    def before_insert(self, mapper, connection, instance):
        instance.unique_id = funcs.hash(instance.sample_id,
                                        instance.sequence)
        instance.junction_num_nts = len(instance.junction_nt)


class Sequence(BaseData):
    """Represents a single unique sequence.

    :param str unique_id: A key over ``sample_id`` and ``sequence`` so the \
        tuple can be maintained unique
    :param str seq_id: A unique identifier for the sequence as output by the \
        sequencer
    :param int sample_id: The ID of the sample from which this sequence came
    :param Relationship sample: Reference to the associated \
        :py:class:`Sample` instance
    :param str alignment: Alignment type (e.g. R1, R2, or R1+R2)
    :param int levenshtein_dist: The optional Levenshtein distance of the \
        sequence to its germline.  Used to identify possible indels and \
        misalignments.
    :param int num_gaps: Number of inserted gaps
    :param int pad_length: The number of pad nucleotides added to the V end \
        of the sequence
    :param int v_match: The number of V-gene nucleotides matching the germline
    :param int v_length: The length of the V-gene segment prior to a streak \
        of mismatches in the CDR3

    :param int j_match: The number of J-gene nucleotides matching the germline
    :param int j_length: The length of the J-gene segment after a streak of \
        mismatches in the CDR3

    :param int pre_cdr3_length: The length of the V-gene prior to the CDR3
    :param int pre_cdr3_match: The number of V-gene nucleotides matching the \
        germline prior to the CDR3

    :param int post_cdr3_length: The length of the J-gene after to the CDR3
    :param int post_cdr3_match: The number of J-gene nucleotides matching the \
        germline after to the CDR3

    :param bool in_frame: If the sequence's CDR3 has a length divisible by 3
    :param bool functional: If the sequence is in-frame and contains no stop \
        codons

    :param bool stop: If the sequence contains a stop codon
    :param int copy_number: Number of reads identical to the sequence in the \
        same sample
    :param str sequence: The (possibly-padded) sequence


    :param str v_gene: The V-gene assigned to the sequence
    :param str j_gene: The J-gene assigned to the sequence
    :param int junction_num_nts: The number of nucleotides in the CDR3
    :param str junction_nt: The nucleotides comprising the CDR3
    :param str junction_aa: The amino-acids comprising the CDR3
    :param str gap_method: The method used to gap the sequence (e.g. IGMT)
    :param str sequence_replaced: The full sequence after being filled in \
        with the germline
    :param str germline: The germline sequence for this sequence

    :param int clone_id: The clone ID to which this sequence belongs
    :param Relationship clone: Reference to the associated :py:class:`Clone` \
        instance

    """
    __tablename__ = 'sequences'
    __table_args__ = (Index('genes', 'v_gene', 'j_gene'),
                      {'mysql_engine': 'TokuDB'})
    __mapper_args__ = {'extension': SequenceExtension()}

    unique_id = Column(String(40), unique=True)

    seq_id = Column(String(128), primary_key=True, index=True)
    sample_id = Column(Integer, ForeignKey(Sample.id), primary_key=True,
                       index=True)
    sample = relationship(Sample, backref=backref('sequences'))

    alignment = Column(String(length=6), index=True)
    probable_indel_or_misalign = Column(Boolean, index=True)

    v_gene = Column(String(512), index=True)
    j_gene = Column(String(512), index=True)

    num_gaps = Column(Integer)
    pad_length = Column(Integer)

    v_match = Column(Integer)
    v_length = Column(Integer)
    j_match = Column(Integer)
    j_length = Column(Integer)

    pre_cdr3_length = Column(Integer)
    pre_cdr3_match = Column(Integer)
    post_cdr3_length = Column(Integer)
    post_cdr3_match = Column(Integer)

    in_frame = Column(Boolean)
    functional = Column(Boolean, index=True)
    stop = Column(Boolean)
    copy_number = Column(Integer, index=True)

    # This is just length(junction_nt) but is included for fast statistics
    # generation over the index
    junction_num_nts = Column(Integer, index=True)

    junction_nt = Column(String(512))
    junction_aa = Column(String(512), index=True)
    gap_method = Column(String(16))

    sequence = Column(String(length=1024), index=True)
    sequence_replaced = Column(String(length=1024), index=True)

    germline = Column(String(length=1024))

    clone_id = Column(Integer, ForeignKey(Clone.id), index=True)
    clone = relationship(Clone, backref=backref('sequences',
                         order_by=seq_id))
    mutations_from_clone = Column(MEDIUMTEXT)


class DuplicateSequence(BaseData):
    """A sequence which is a duplicate of a :py:class:`Sequence`.  This is
    used to minimize the size of the sequences table.  The ``copy_number``
    attribute of :py:class:`Sequence` instances is equal to the number of
    its duplicate sequences plus one.

    :param str duplicate_seq_id: The identifier of the sequence in the same \
        sample with the same sequence
    :param Relationship duplicate_seq: Reference to the associated \
        :py:class:`Sequence` instance of which this is a duplicate

    :param int sample_id: The ID of the sample from which this sequence came
    :param Relationship sample: Reference to the associated \
        :py:class:`Sample` instance

    :param str seq_id: A unique identifier for the sequence as output by the \
        sequencer

    """
    __tablename__ = 'duplicate_sequences'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    seq_id = Column(String(length=128), primary_key=True)

    duplicate_seq_id = Column(String(length=128),
                              ForeignKey('sequences.seq_id'),
                              primary_key=True,
                              index=True)
    duplicate_seq = relationship(Sequence,
                                 backref=backref('duplicate_sequences',
                                                 order_by=duplicate_seq_id))

    sample_id = Column(Integer, ForeignKey(Sample.id),
                       primary_key=True)
    sample = relationship(Sample, backref=backref('duplicate_sequences',
                          order_by=seq_id))


class NoResult(BaseData):
    """A sequence which could not be match with a V or J.

    :param str seq_id: A unique identifier for the sequence as output by the \
        sequencer
    :param int sample_id: The ID of the sample from which this sequence came
    :param Relationship sample: Reference to the associated \
        :py:class:`Sample` instance
    :param str sequence: The sequence of the non-identifiable input

    """
    __tablename__ = 'noresults'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    seq_id = Column(String(length=128), primary_key=True)

    sample_id = Column(Integer, ForeignKey(Sample.id),
                       primary_key=True)
    sample = relationship(Sample, backref=backref('noresults',
                          order_by=seq_id))

    sequence = Column(String(length=1024))


class ModificationLog(BaseData):
    __tablename__ = 'modification_logs'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    id = Column(Integer, primary_key=True)
    datetime = Column(DateTime, default=datetime.datetime.utcnow)

    action_type = Column(String(length=128))
    info = Column(String(length=1024))
