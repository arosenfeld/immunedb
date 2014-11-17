from sqlalchemy import Column, Boolean, Integer, String, Text, Date, \
    ForeignKey, UniqueConstraint, Index, func
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.dialects.mysql import TEXT, MEDIUMTEXT

from sldb.common.settings import DATABASE_SETTINGS

assert DATABASE_SETTINGS['master_metadata'] is not None, ('Must specify master '
    'metadata')
assert DATABASE_SETTINGS['data_metadata'] is not None, ('Must specify data '
    'metadata')

BaseMaster = declarative_base(
    metadata=DATABASE_SETTINGS['master_metadata'])
BaseData = declarative_base(
    metadata=DATABASE_SETTINGS['data_metadata'])

class Study(BaseMaster):
    """Represents a high-level study (e.g. Lupus)"""
    __tablename__ = 'studies'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    id = Column(Integer, primary_key=True)
    # The name of the study
    name = Column(String(length=128), unique=True)
    # Some arbitrary information if necessary
    info = Column(String(length=1024))


class Subject(BaseMaster):
    __tablename__ = 'subjects'
    __table_args__ = (UniqueConstraint('study_id', 'identifier'),
                      {'mysql_engine': 'TokuDB'})

    id = Column(Integer, primary_key=True)

    identifier = Column(String(64))
    study_id = Column(Integer, ForeignKey(Study.id))
    study = relationship(Study, backref=backref('subjects',
                         order_by=identifier))


class Sample(BaseMaster):
    """A sample from a study.  Generally this is from a single subject and
    tissue"""
    __tablename__ = 'samples'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    # Unique ID for the sample
    id = Column(Integer, primary_key=True)
    # A name for the sample
    name = Column(String(128), unique=True)
    # Some arbitrary information if necessary
    info = Column(String(1024))

    # The date the sample was taken
    date = Column(Date)

    # Reference to the study
    study_id = Column(Integer, ForeignKey(Study.id))
    study = relationship(Study, backref=backref('samples', order_by=(date,
                                                name)))

    subject_id = Column(Integer, ForeignKey(Subject.id), index=True)
    subject = relationship(Subject, backref=backref('samples',
                           order_by=(id)))
    subset = Column(String(16))
    tissue = Column(String(16))
    disease = Column(String(32))
    lab = Column(String(128))
    experimenter = Column(String(128))

    # Number of valid sequences in the sample
    valid_cnt = Column(Integer, default=0)
    # Number of invalid sequences in the sample
    no_result_cnt = Column(Integer, default=0)
    # Number of functional sequences in the sample
    functional_cnt = Column(Integer, default=0)


class SampleStats(BaseData):
    """Aggregate statistics for a sample.  This exists to reduce the time
    queries take for a sample."""
    __tablename__ = 'sample_stats'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    # Reference to the sample the stats are for
    sample_id = Column(Integer, ForeignKey(Sample.id),
                       primary_key=True)
    sample = relationship(Sample, backref=backref('sample_stats',
                          order_by=sample_id))

    # The filter type for the stats (e.g. unique, clones)
    filter_type = Column(String(length=255), primary_key=True)

    # Distributions stored as JSON for a given field in the sample
    v_match_dist = Column(MEDIUMTEXT)
    v_length_dist = Column(MEDIUMTEXT)

    j_match_dist = Column(MEDIUMTEXT)
    j_length_dist = Column(MEDIUMTEXT)

    v_call_dist = Column(MEDIUMTEXT)
    j_call_dist = Column(MEDIUMTEXT)

    v_gap_length_dist = Column(MEDIUMTEXT)
    j_gap_length_dist = Column(MEDIUMTEXT)

    copy_number_iden_dist = Column(MEDIUMTEXT)
    cdr3_length_dist = Column(MEDIUMTEXT)

    # Total sequences in filter
    sequence_cnt = Column(Integer)
    # Number of in-frame sequences
    in_frame_cnt = Column(Integer)
    # Number of sequences with a stop codon
    stop_cnt = Column(Integer)


class CloneGroup(BaseMaster):
    """A clone which is dictated by V, J, CDR3."""
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
    """A group of sequences likely originating from the same germline"""
    __tablename__ = 'clones'
    __table_args__ = (Index('v_gene', 'j_gene', 'cdr3_num_nts', 'subject_id'),
                      {'mysql_engine': 'TokuDB'})

    id = Column(Integer, primary_key=True)
    cdr3_nt = Column(String(length=512))

    # These are necessary during creation of clones, but will
    # always be redundant with the associated grouping after completion
    v_gene = Column(String(length=512))
    j_gene = Column(String(length=128))
    cdr3_num_nts = Column(Integer)
    subject_id = Column(Integer, ForeignKey(Subject.id))
    #

    group_id = Column(Integer, ForeignKey(CloneGroup.id),
                      index=True)
    group = relationship(CloneGroup, backref=backref('clones',
                         order_by=(id)))


class Sequence(BaseData):
    """Represents a single unique sequence."""
    __tablename__ = 'sequences'
    __table_args__ = (Index('v_call', 'j_call',),
                      {'mysql_engine': 'TokuDB'})

    seq_id = Column(String(length=128), primary_key=True)

    v_call = Column(String(512))
    j_call = Column(String(512))

    junction_nt = Column(String(512))
    junction_aa = Column(String(512))
    gap_method = Column(String(16))

    sequence_replaced = Column(String(length=1024), unique=True, index=True)

    clone_id = Column(Integer, ForeignKey(Clone.id), index=True)
    clone = relationship(Clone, backref=backref('sequences',
                         order_by=seq_id))


class SequenceMapping(BaseData):
    __tablename__ = 'sequence_mapping'
    __table_args__ = (UniqueConstraint('seq_id', 'sample_id'),
                      {'mysql_engine': 'TokuDB'})

    identity_seq_id = Column(String(128), ForeignKey(Sequence.seq_id),
                             primary_key=True)
    identity_seq = relationship(Sequence, backref=backref('mapping'))

    seq_id = Column(String(128))
    sample_id = Column(Integer, ForeignKey(Sample.id), index=True)
    sample = relationship(Sample, backref=backref('mapping'))

    alignment = Column(String(length=6))
    copy_number = Column(Integer)

    v_match = Column(Integer)
    v_length = Column(Integer)
    j_match = Column(Integer)
    j_length = Column(Integer)

    in_frame = Column(Boolean)
    functional = Column(Boolean, index=True)
    stop = Column(Boolean)

    sequence = Column(String(length=1024))


class CloneFrequency(BaseData):
    """Frequency statistics for a clone with different filters."""
    __tablename__ = 'clone_frequencies'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    sample_id = Column(Integer, ForeignKey(Sample.id),
                       primary_key=True)
    sample = relationship(Sample, backref=backref('clone_frequencies',
                          order_by=sample_id))

    clone_id = Column(Integer, ForeignKey(Clone.id),
                        primary_key=True)
    clone = relationship(Clone, backref=backref('clone_frequencies',
                         order_by=sample_id))

    filter_type = Column(String(length=255), primary_key=True)

    unique_sequences = Column(Integer)
    total_sequences = Column(Integer)


class NoResult(BaseData):
    """A sequence which could not be match with a V or J."""
    __tablename__ = 'noresults'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    seq_id = Column(String(length=128), primary_key=True)

    sample_id = Column(Integer, ForeignKey(Sample.id),
                       primary_key=True)
    sample = relationship(Sample, backref=backref('noresults',
                          order_by=seq_id))
