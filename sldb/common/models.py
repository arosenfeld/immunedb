from sqlalchemy import Column, Boolean, Integer, String, Text, Date, \
    ForeignKey, UniqueConstraint, Index
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref
from sqlalchemy.ext.declarative import declared_attr
from sqlalchemy.dialects.mysql import TEXT, MEDIUMTEXT


Base = declarative_base()


class Study(Base):
    __tablename__ = 'studies'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    id = Column(Integer, primary_key=True)
    name = Column(String(length=128), unique=True)
    info = Column(String(length=1024))


class Sample(Base):
    __tablename__ = 'samples'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    id = Column(Integer, primary_key=True)
    name = Column(String(128))
    info = Column(String(1024))

    date = Column(Date)

    study_id = Column(Integer, ForeignKey('studies.id'))
    study = relationship('Study', backref=backref('samples', order_by=(date,
                                                  name)))

    valid_cnt = Column(Integer)
    no_result_cnt = Column(Integer)
    functional_cnt = Column(Integer)


class SampleStats(Base):
    __tablename__ = 'sample_stats'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    sample_id = Column(Integer, ForeignKey('samples.id'),
                       primary_key=True)
    sample = relationship('Sample', backref=backref('sample_stats',
                          order_by=sample_id))

    filter_type = Column(String(length=255), primary_key=True)

    sequence_cnt = Column(Integer)
    v_match_dist = Column(MEDIUMTEXT)
    v_length_dist = Column(MEDIUMTEXT)

    j_match_dist = Column(MEDIUMTEXT)
    j_length_dist = Column(MEDIUMTEXT)

    v_call_dist = Column(MEDIUMTEXT)
    j_call_dist = Column(MEDIUMTEXT)

    v_gap_length_dist = Column(MEDIUMTEXT)
    j_gap_length_dist = Column(MEDIUMTEXT)

    copy_number_close_dist = Column(MEDIUMTEXT)
    copy_number_iden_dist = Column(MEDIUMTEXT)
    cdr3_len_dist = Column(MEDIUMTEXT)

    clone_dist = Column(MEDIUMTEXT)

    in_frame_cnt = Column(Integer)
    stop_cnt = Column(Integer)
    mutation_inv_cnt = Column(Integer)


class Sequence(Base):
    __tablename__ = 'sequences'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    seq_id = Column(String(length=128), primary_key=True)

    sample_id = Column(Integer, ForeignKey('samples.id'),
                       primary_key=True)
    sample = relationship('Sample', backref=backref('sequences',
                          order_by=seq_id))
    order = Column(Integer)

    functional = Column(Boolean)
    in_frame = Column(Boolean)
    stop = Column(Boolean)
    mutation_invariate = Column(Boolean)  # TODO: Is this a typo?

    v_match = Column(Integer)
    v_length = Column(Integer)

    j_match = Column(Integer)
    j_length = Column(Integer)

    v_call = Column(String(512))
    j_call = Column(String(512))

    v_gap_length = Column(Integer)
    j_gap_length = Column(Integer)
    # No junction_gap_length; use len(junction_nt)
    junction_nt = Column(String(512))
    junction_aa = Column(String(512))
    gap_method = Column(String(16))

    subject = Column(String(128))
    subset = Column(String(16))

    tissue = Column(String(16))

    disease = Column(String(32))

    date = Column(Date)

    lab = Column(String(128))
    experimenter = Column(String(128))

    copy_number_close = Column(Integer)
    collapse_to_close = Column(Integer)
    copy_number_iden = Column(Integer)
    collapse_to_iden = Column(Integer)

    sequence = Column(String(length=1024))
    sequence_replaced = Column(String(length=1024))

    clone_id = Column(Integer, ForeignKey('clones.id'),
                      index=True)
    clone = relationship('Clone', backref=backref('sequences',
                         order_by=seq_id))

    clone_size = Column(Integer)
    clone_copy_number = Column(Integer)


class Clone(Base):
    __tablename__ = 'clones'
    __table_args__ = (UniqueConstraint('v_gene', 'j_gene', 'cdr3',
                      'cdr3_num_nts'),
                      {'mysql_engine': 'TokuDB'})

    id = Column(Integer, primary_key=True)

    v_gene = Column(String(length=200))
    j_gene = Column(String(length=200))
    cdr3 = Column(String(length=128))
    cdr3_num_nts = Column(Integer)

    germline = Column(String(length=1024))


class CloneFrequency(Base):
    __tablename__ = 'clone_frequencies'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    sample_id = Column(Integer, ForeignKey('samples.id'),
                       primary_key=True)
    sample = relationship('Sample', backref=backref('clone_frequencies',
                          order_by=sample_id))

    clone_id = Column(Integer, ForeignKey('clones.id'),
                      primary_key=True)
    clone = relationship('Clone', backref=backref('clone_frequencies',
                                                  order_by=sample_id))

    filter_type = Column(String(length=255), primary_key=True)

    size = Column(Integer)
    copy_number = Column(Integer)


class NoResult(Base):
    __tablename__ = 'noresults'
    __table_args__ = {'mysql_engine': 'TokuDB'}

    seq_id = Column(String(length=128), primary_key=True)

    sample_id = Column(Integer, ForeignKey('samples.id'),
                       primary_key=True)
    sample = relationship('Sample', backref=backref('noresults',
                          order_by=seq_id))
