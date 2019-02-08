import datetime
import json

from sqlalchemy import (Column, Boolean, Float, Integer, String, DateTime,
                        ForeignKey, UniqueConstraint, Index, event)
from sqlalchemy.dialects.mysql import MEDIUMTEXT
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import ColumnProperty, relationship, backref
from sqlalchemy.orm.session import Session
from sqlalchemy.schema import ForeignKeyConstraint, PrimaryKeyConstraint

import immunedb.util.funcs as funcs

Base = declarative_base()
MAX_CDR3_NTS = 96
MAX_CDR3_AAS = int(MAX_CDR3_NTS / 3)
MAX_SEQ_LEN = 512
MAX_GENE_LEN = 64
MAX_INDEL_LEN = 64
CDR3_OFFSET = 309


def deserialize_gaps(gaps):
    if gaps is None:
        return []
    return [
        [int(p) for p in g.split('-')]
        for g in gaps.split(',')
    ]


def serialize_gaps(gaps):
    if gaps is None or len(gaps) == 0:
        return None
    return ','.join(
        ['{}-{}'.format(start, end) for (start, end) in sorted(gaps)]
    )


class Study(Base):
    """A study which aggregates related samples.

    :param int id: An auto-assigned unique identifier for the study
    :param str name: A unique name for the study
    :param str info: Optional information about the study

    """
    __tablename__ = 'studies'
    __table_args__ = {'mysql_row_format': 'DYNAMIC'}

    id = Column(Integer, primary_key=True)
    name = Column(String(length=128), unique=True)
    info = Column(String(length=1024))


class Subject(Base):
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
                      {'mysql_row_format': 'DYNAMIC'})

    id = Column(Integer, primary_key=True)

    identifier = Column(String(64))
    study_id = Column(Integer, ForeignKey(Study.id))
    study = relationship(Study, backref=backref('subjects',
                         order_by=identifier))


class Sample(Base):
    """A sample of sequences.

    :param int id: An auto-assigned unique identifier for the sample
    :param str name: A unique name for the sample as defined by the \
        experimenter

    :param int study_id: The ID of the study under which the subject was \
        sampled
    :param Relationship study: Reference to the associated :py:class:`Study` \
        instance

    :param int subject_id: The ID of the subject from which the sample was \
        taken
    :param Relationship subject: Reference to the associated \
        :py:class:`Subject` instance


    :param float v_ties_mutations: Average mutation rate of sequences in the \
        sample
    :param float v_ties_len: Average length of sequences in the sample

    """
    __tablename__ = 'samples'
    __table_args__ = {'mysql_row_format': 'DYNAMIC'}

    id = Column(Integer, primary_key=True)
    name = Column(String(128), unique=True)

    study_id = Column(Integer, ForeignKey(Study.id))
    study = relationship(Study, backref=backref('samples', order_by=name))

    subject_id = Column(Integer, ForeignKey(Subject.id), index=True)
    subject = relationship(Subject, backref=backref('samples',
                           order_by=(id)))

    v_ties_mutations = Column(Float)
    v_ties_len = Column(Float)

    @property
    def metadata_dict(self):
        return {m.key: m.value for m in self.metadata_models}

    @property
    def stats(self):
        for s in self.sample_stats:
            if s.filter_type == 'all' and s.outliers and not s.full_reads:
                return s


class SampleMetadata(Base):
    __tablename__ = 'sample_metadata'
    __table_args__ = {'mysql_row_format': 'DYNAMIC'}

    sample_id = Column(Integer, ForeignKey(Sample.id, ondelete='CASCADE'),
                       primary_key=True, nullable=False)
    sample = relationship(Sample, backref=backref('metadata_models',
                          order_by=sample_id))

    key = Column(String(length=32), primary_key=True, index=True,
                 nullable=False)
    value = Column(String(length=64), index=True)


class SampleStats(Base):
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

    :param str v_identity_dist: Distribution of V gene identity
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
    __table_args__ = (
        Index('stat_cover', 'sample_id', 'outliers', 'full_reads',
              'filter_type', 'sequence_cnt', 'in_frame_cnt', 'stop_cnt',
              'functional_cnt', 'no_result_cnt'),
        {'mysql_row_format': 'DYNAMIC'})

    sample_id = Column(Integer, ForeignKey(Sample.id),
                       primary_key=True)
    sample = relationship(Sample, backref=backref('sample_stats',
                          order_by=sample_id, cascade='all, delete-orphan'))

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

    sp_fwr_dist = Column(MEDIUMTEXT)
    sp_cdr_dist = Column(MEDIUMTEXT)

    quality_dist = Column(MEDIUMTEXT)

    sequence_cnt = Column(Integer)
    in_frame_cnt = Column(Integer)
    stop_cnt = Column(Integer)
    functional_cnt = Column(Integer)
    no_result_cnt = Column(Integer)


class Clone(Base):
    """A group of sequences likely originating from the same germline

    :param int id: An auto-assigned unique identifier for the clone
    :param bool functional: If the clone is functional
    :param str v_gene: The V-gene assigned to the sequence
    :param str j_gene: The J-gene assigned to the sequence
    :param str cdr3_nt: The consensus nucleotides for the clone
    :param int cdr3_num_nts: The number of nucleotides in the group's CDR3
    :param str cdr3_aa: The amino-acid sequence of the group's CDR3
    :param int subject_id: The ID of the subject to which the clone belongs
    :param Relationship subject: Reference to the associated \
        :py:class:`Subject` instance
    :param str germline: The germline sequence for this clone
    :param str tree: The textual representation of the clone's lineage tree
    :param int parent_id: The (possibly null) ID of the clone's parent

    """
    __tablename__ = 'clones'
    __table_args__ = (Index('bucket', 'v_gene', 'j_gene', 'subject_id',
                            'cdr3_num_nts', 'insertions', 'deletions'),
                      Index('aa_bucket', 'v_gene', 'j_gene', 'subject_id',
                            'cdr3_aa'),
                      Index('tree', 'tree', mysql_length=1),
                      {'mysql_row_format': 'DYNAMIC'})
    id = Column(Integer, primary_key=True)

    functional = Column(Boolean, index=True)

    v_gene = Column(String(length=MAX_GENE_LEN), index=True)
    j_gene = Column(String(length=MAX_GENE_LEN), index=True)

    _insertions = Column('insertions', String(MAX_INDEL_LEN), index=True)
    _deletions = Column('deletions', String(MAX_INDEL_LEN), index=True)

    cdr3_nt = Column(String(length=MAX_CDR3_NTS))
    cdr3_num_nts = Column(Integer, index=True)
    cdr3_aa = Column(String(length=MAX_CDR3_AAS))

    subject_id = Column(Integer, ForeignKey(Subject.id), index=True)
    subject = relationship(Subject, backref=backref('clones',
                           order_by=(v_gene, j_gene, cdr3_num_nts, cdr3_aa)))

    germline = Column(String(length=MAX_SEQ_LEN))
    tree = Column(MEDIUMTEXT)

    overall_unique_cnt = Column(Integer, index=True)  # Denormalized
    overall_instance_cnt = Column(Integer, index=True)  # Denormalized
    overall_total_cnt = Column(Integer, index=True)  # Denormalized

    parent_id = Column(Integer, ForeignKey('clones.id'), index=True)

    children = relationship('Clone')
    parent = relationship('Clone', remote_side=[id], back_populates='children')

    @hybrid_property
    def insertions(self):
        """Returns the list of insertion position/length pairs"""
        return deserialize_gaps(self._insertions)

    @insertions.setter
    def insertions(self, value):
        """Sets the list of insertion position/length pairs"""
        self._insertions = serialize_gaps(value)

    @hybrid_property
    def deletions(self):
        """Returns the list of deletion position/length pairs"""
        return deserialize_gaps(self._deletions)

    @deletions.setter
    def deletions(self, value):
        """Sets the list of deletions position/length pairs"""
        self._deletions = serialize_gaps(value)

    @property
    def regions(self):
        """Returns the IMGT region boundaries for the clone"""
        regions = funcs.get_regions(self.insertions)
        regions.append(self.cdr3_num_nts)
        regions.append(len(self.germline) - sum(regions))
        return regions

    @property
    def consensus_germline(self):
        """Returns the consensus germline for the clone"""
        cdr3_start = CDR3_OFFSET
        if self.insertions is not None:
            cdr3_start += sum((e[1] for e in self.insertions))
        return ''.join([
            self.germline[0:cdr3_start],
            self.cdr3_nt,
            self.germline[cdr3_start + self.cdr3_num_nts:]
        ])

    @property
    def overall_unique_cnt_with_subclones(self):
        return self.overall_unique_cnt + sum([
            s.overall_unique_cnt for s in self.children])

    @property
    def overall_total_cnt_with_subclones(self):
        return self.overall_total_cnt + sum([
            s.overall_total_cnt for s in self.children])

    @property
    def overall_instance_cnt_with_subclones(self):
        return self.overall_instance_cnt + sum([
            s.overall_instance_cnt for s in self.children])

    @property
    def overall_stats(self):
        return [s for s in self.stats if not s.sample_id][0]


class CloneStats(Base):
    """Stores statistics for a given clone and sample.  If sample is null the
    statistics are for the specified clone in all samples.

    :param int clone_id: The clone ID
    :param Relationship clone: Reference to the associated \
        :py:class:`Clone` instance

    :param bool functional: If the associated clone is functional.  This is a \
        denormalized field.

    :param int sample_id: The sample ID
    :param Relationship sample: Reference to the associated \
        :py:class:`Sample` instance

    :param int unique_cnt: The number of unique sequences in the clone in the \
        sample
    :param int total_cnt: The number of total sequences in the clone in the \
        sample

    :param str mutations: A JSON stanza of mutation count information

    """
    __tablename__ = 'clone_stats'
    __table_args__ = (
        Index('stats_cover', 'sample_id', 'clone_id', 'functional',
              'total_cnt', 'unique_cnt'),
        {'mysql_row_format': 'DYNAMIC'})

    id = Column(Integer, primary_key=True)
    clone_id = Column(Integer, ForeignKey(Clone.id, ondelete='CASCADE'))
    clone = relationship(Clone, backref=backref('stats'))
    subject_id = Column(
        Integer, ForeignKey(Subject.id, ondelete='CASCADE'))  # Denormalized

    functional = Column(Boolean, index=True)  # Denormalized

    sample_id = Column(Integer, ForeignKey(Sample.id, ondelete='CASCADE'))
    sample = relationship(Sample, backref=backref('clone_stats'))

    unique_cnt = Column(Integer)
    total_cnt = Column(Integer)

    mutations = Column(MEDIUMTEXT)
    avg_v_identity = Column(Float)

    top_copy_seq_ai = Column(Integer, ForeignKey('sequences.ai'))
    top_copy_seq_sequence = Column(String(length=MAX_SEQ_LEN))
    top_copy_seq_copies = Column(Integer)
    top_copy_seq = relationship('Sequence')

    @property
    def v_mutations(self):
        muts = json.loads(self.mutations).get('regions', {})
        aggregate = 0
        for region in ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3']:
            region_muts = muts.get(region, {})
            aggregate += sum((
                sum([s['total'] for s in type_muts])
                for type_muts in region_muts.values()
            ))
        return aggregate


class SelectionPressure(Base):
    __tablename__ = 'selection_pressure'
    __table_args__ = (
        UniqueConstraint('clone_id', 'sample_id', 'threshold'),
        {'mysql_row_format': 'DYNAMIC'})
    id = Column(Integer, primary_key=True)
    clone_id = Column(Integer, ForeignKey(Clone.id, ondelete='CASCADE'))
    clone = relationship(Clone, backref=backref('selection_pressure'))

    sample_id = Column(Integer, ForeignKey(Sample.id))
    sample = relationship(Sample, backref=backref('selection_pressure'))

    threshold = Column(String(length=10), index=True)

    expected_fwr_s = Column(Float)
    expected_cdr_s = Column(Float)
    expected_fwr_r = Column(Float)
    expected_cdr_r = Column(Float)

    observed_fwr_s = Column(Float)
    observed_cdr_s = Column(Float)
    observed_fwr_r = Column(Float)
    observed_cdr_r = Column(Float)

    sigma_fwr = Column(Float)
    sigma_cdr = Column(Float)

    sigma_fwr_cilower = Column(Float)
    sigma_fwr_ciupper = Column(Float)
    sigma_cdr_cilower = Column(Float)
    sigma_cdr_ciupper = Column(Float)

    sigma_p_fwr = Column(Float)
    sigma_p_cdr = Column(Float)

    def to_dict(self):
        fields = [
            'expected_fwr_s',
            'expected_cdr_s',
            'expected_fwr_r',
            'expected_cdr_r',

            'observed_fwr_s',
            'observed_cdr_s',
            'observed_fwr_r',
            'observed_cdr_r',

            'sigma_fwr',
            'sigma_cdr',

            'sigma_fwr_cilower',
            'sigma_fwr_ciupper',
            'sigma_cdr_cilower',
            'sigma_cdr_ciupper',

            'sigma_p_fwr',
            'sigma_p_cdr'
        ]
        return {field: getattr(self, field) for field in fields}


class Sequence(Base):
    """Represents a single unique sequence.

    :param int ai: An auto-incremented value for the sequence

    :param int subject_id: The ID of the subject for this subject
    :param str seq_id: A unique identifier for the sequence as output by the \
        sequencer

    :param int sample_id: The ID of the sample from which this sequence came
    :param Relationship sample: Reference to the associated \
        :py:class:`Sample` instance

    :param bool partial: If the sequence is a partial read
    :param bool probable_indel_or_misalign: If the sequence likely has an \
        indel or is a bad alignment

    :param str v_gene: The V-gene assigned to the sequence
    :param str j_gene: The J-gene assigned to the sequence

    :param int num_gaps: Number of inserted gaps within the V read
    :param int seq_start: The offset from the germline where the sequence \
        starts

    :param int v_match: The number of V-gene nucleotides matching the germline
    :param int v_length: The length of the V-gene segment prior to a streak \
        of mismatches in the CDR3
    :param int j_match: The number of J-gene nucleotides matching the germline
    :param int j_length: The length of the J-gene segment after a streak of \
        mismatches in the CDR3

    :param str removed_prefix: The sequence (if any) which was removed from \
        the beginning of the sequence during alignment.  Possibly used \
        during indel correction
    :param str removed_prefix_qual: The quality (if any) which was removed \
        from the beginning of the sequence during alignment.  Possibly used \
        during indel correction

    :param int pre_cdr3_length: The length of the V-gene prior to the CDR3
    :param int pre_cdr3_match: The number of V-gene nucleotides matching the \
        germline prior to the CDR3

    :param int post_cdr3_length: The length of the J-gene after to the CDR3
    :param int post_cdr3_match: The number of J-gene nucleotides matching the \
        germline after to the CDR3

    :param bool in_frame: If the sequence's CDR3 has a length divisible by 3
    :param bool functional: If the sequence is functional
    :param bool stop: If the sequence contains a stop codon
    :param int copy_number: Number of reads in the sample which collapsed to \
        this sequence

    :param int cdr3_num_nts: The number of nucleotides in the CDR3
    :param str cdr3_nt: The nucleotides comprising the CDR3
    :param str cdr3_aa: The amino-acids comprising the CDR3

    :param str sequence: The (possibly-padded) sequence
    :param str quality: Optional Phred quality score (in Sanger format) for \
        each base in ``sequence``

    :param str germline: The germline sequence for this sequence

    :param int clone_id: The clone ID to which this sequence belongs
    :param Relationship clone: Reference to the associated :py:class:`Clone` \
        instance
    :param str mutations_from_clone: A JSON stanza with mutation information


    """
    __tablename__ = 'sequences'
    __table_args__ = (
        Index('local_align_bucket', 'sample_id', 'locally_aligned', 'v_gene',
              'j_gene', 'cdr3_num_nts', 'copy_number', 'ai'),
        Index('subject_bucket', 'subject_id', 'v_gene', 'j_gene',
              'cdr3_num_nts', 'insertions', 'deletions'),
        Index('subject_clone_bucket', 'subject_id', 'clone_id', 'v_gene',
              'j_gene', 'cdr3_num_nts', 'insertions', 'deletions'),
        Index('locally_aligned', 'locally_aligned', 'v_gene', 'j_gene',
              'cdr3_num_nts', 'sample_id'),
        UniqueConstraint('sample_id', 'seq_id'),
        PrimaryKeyConstraint('sample_id', 'ai'),
        # NOTE: This is because sqlalchemy doesn't properly preserve the column
        # ordering.
        UniqueConstraint('sample_id', 'ai'),
        {'mysql_row_format': 'DYNAMIC'}
    )

    def __init__(self, **kwargs):
        self.insertions = kwargs.pop('insertions', None)
        self.deletions = kwargs.pop('deletions', None)
        super(Sequence, self).__init__(**kwargs)

    sample_id = Column(Integer, ForeignKey(Sample.id))
    ai = Column(Integer, autoincrement=True, unique=True)

    subject_id = Column(Integer, ForeignKey(Subject.id), index=True)
    subject = relationship(Subject)

    seq_id = Column(String(64), index=True)
    sample = relationship(Sample, backref=backref('sequences'))

    partial = Column(Boolean, index=True)
    rev_comp = Column(Boolean)

    probable_indel_or_misalign = Column(Boolean)
    locally_aligned = Column(Boolean, default=False, nullable=False)

    _deletions = Column('deletions', String(MAX_INDEL_LEN))
    _insertions = Column('insertions', String(MAX_INDEL_LEN))

    v_gene = Column(String(MAX_GENE_LEN))
    j_gene = Column(String(MAX_GENE_LEN))

    num_gaps = Column(Integer)
    seq_start = Column(Integer)

    v_match = Column(Integer)
    v_length = Column(Integer)
    j_match = Column(Integer)
    j_length = Column(Integer)

    removed_prefix = Column(String(256))
    removed_prefix_qual = Column(String(256))
    v_mutation_fraction = Column(Float)

    pre_cdr3_length = Column(Integer)
    pre_cdr3_match = Column(Integer)
    post_cdr3_length = Column(Integer)
    post_cdr3_match = Column(Integer)

    in_frame = Column(Boolean)
    functional = Column(Boolean)
    stop = Column(Boolean)
    copy_number = Column(Integer, server_default='0', nullable=False,
                         index=True)

    # This is just length(cdr3_nt) but is included for fast statistics
    # generation over the index
    cdr3_num_nts = Column(Integer)

    cdr3_nt = Column(String(MAX_CDR3_NTS))
    cdr3_aa = Column(String(MAX_CDR3_AAS))

    sequence = Column(String(length=MAX_SEQ_LEN))
    quality = Column(String(length=MAX_SEQ_LEN))

    germline = Column(String(length=MAX_SEQ_LEN))

    clone_id = Column(Integer, ForeignKey(Clone.id, ondelete='SET NULL'),
                      index=True)
    clone = relationship(Clone, backref=backref('sequences',
                         order_by=seq_id))
    mutations_from_clone = Column(MEDIUMTEXT)

    @hybrid_property
    def deletions(self):
        """Returns the list of deletion position/length pairs"""
        return deserialize_gaps(self._deletions)

    @deletions.setter
    def deletions(self, value):
        """Sets the list of deletions position/length pairs"""
        self._deletions = serialize_gaps(value)

    @hybrid_property
    def insertions(self):
        """Returns the list of insertions position/length pairs"""
        return deserialize_gaps(self._insertions)

    @insertions.setter
    def insertions(self, value):
        """Sets the list of insertions position/length pairs"""
        self._insertions = serialize_gaps(value)

    @property
    def original_sequence(self):
        """Returns the original sequence given with the J end trimmed to the
           germline

        """
        return '{}{}'.format(
            self.removed_prefix,
            self.sequence.replace('-', '')
        )

    @property
    def original_quality(self):
        """Returns the original quality given with the J end trimmed to the
           germline

        """
        if self.quality is None:
            return None
        return '{}{}'.format(self.removed_prefix_qual or '',
                             self.quality.replace(' ', ''))

    @property
    def clone_sequence(self):
        """Gets the sequence within the context of the associated clone by
        adding insertions from other sequences to this one.

        """
        seq = self.sequence
        if self.clone is None:
            return seq
        for ins in sorted(self.clone.insertions):
            if ins in self.insertions:
                continue
            pos, size = ins
            seq = seq[:pos] + ('-' * size) + seq[pos:]

        return seq

    @property
    def regions(self):
        """Returns the IMGT region boundaries for the sequence"""
        regions = funcs.get_regions(self.insertions)
        regions.append(self.cdr3_num_nts)
        regions.append(len(self.germline) - sum(regions))
        return regions

    @property
    def v_cigar(self):
        cigar = []
        if self.removed_prefix:
            cigar.extend((len(self.removed_prefix), 'S'))
        if self.seq_start:
            cigar.extend((self.seq_start, 'N'))
        end = self.seq_start + self.pre_cdr3_length + self.num_gaps
        ref = self.germline[self.seq_start:end]
        qry = self.sequence[self.seq_start:end]

        cigar.extend(funcs.get_cigar(ref, qry))
        return ''.join([str(s) for s in cigar])

    @property
    def j_cigar(self):
        return ''.join(funcs.get_cigar(
            self.germline[-self.post_cdr3_length:],
            self.sequence[-self.post_cdr3_length:]
        ))

    @property
    def germline_d_masked(self):
        v_end = self.seq_start + self.pre_cdr3_length + self.num_gaps
        return ''.join((
            self.germline[:v_end],
            'N' * self.cdr3_num_nts,
            self.germline[-self.post_cdr3_length:]
        ))

    def get_v_extent(self, in_clone):
        """Returns the estimated V length, including the portion in the
           CDR3

        """
        extent = self.v_length + self.num_gaps + self.seq_start

        if not in_clone:
            return extent

        return extent + len(self.clone_sequence) - len(self.sequence)


class NoResult(Base):
    """A sequence which could not be match with a V or J.

    :param int pk: A primary key for this no result
    :param str seq_id: A unique identifier for the sequence as output by the \
        sequencer
    :param int sample_id: The ID of the sample from which this sequence came
    :param Relationship sample: Reference to the associated \
        :py:class:`Sample` instance
    :param str sequence: The sequence of the non-identifiable input
    :param str sequence: The quality of the non-identifiable input

    """
    __tablename__ = 'noresults'
    __table_args__ = (
        Index('sample_seq_id', 'sample_id', 'seq_id'),
        {'mysql_row_format': 'DYNAMIC'}
    )

    pk = Column(Integer, primary_key=True)

    seq_id = Column(String(length=64))

    sample_id = Column(Integer, ForeignKey(Sample.id, ondelete='CASCADE'))
    sample = relationship(Sample)

    # Allow longer sequences here since they aren't aligned and we don't know
    # the length.
    sequence = Column(String(length=MAX_SEQ_LEN * 2))
    quality = Column(String(length=MAX_SEQ_LEN * 2))

    reason = Column(String(256))


class ModificationLog(Base):
    """A log message for a database modification

    :param int id: The ID of the log message
    :param datetime datetime: The date and time of the message
    :param str action_type: A short string representing the action
    :param str info: A JSON stanza with log message information

    """
    __tablename__ = 'modification_logs'
    __table_args__ = {'mysql_row_format': 'DYNAMIC'}

    id = Column(Integer, primary_key=True)
    datetime = Column(DateTime, default=datetime.datetime.utcnow)

    action_type = Column(String(length=128))
    info = Column(String(length=1024))


class SequenceCollapse(Base):
    """A one to many table that links sequence from different samples that
    collapse to one another.  This is used instead of a field in
    :py:class:`Sequence` for performance reasons.

    :param int sample_id: The ID of the sample with the sequence being
       collapsed
    :param int seq_ai: The auto-increment value of the sequence being
       collapsed
    :param Relationship clone: Reference to the associated \
        :py:class:`Sequence` instance being collapsed

    :param int collapse_to_subject_sample_id: The ID of the sample in which \
        the collapse to sequence belongs
    :param int collapse_to_subject_seq_ai: The auto-increment value of the \
        sequence collapsing to
    :param int collapse_to_subject_seq_id: The sequence ID of the \
        sequence collapsing to.  This is a denormalized field.
    :param int instances_in_subject: The number of instance of the sequence \
        in the subject
    :param int copy_number_in_subject: The aggregate copy number of the \
        sequence in the subject

    """
    __tablename__ = 'sequence_collapse'
    __table_args__ = (
        PrimaryKeyConstraint('sample_id', 'seq_ai'),
        ForeignKeyConstraint(
            ['sample_id', 'seq_ai'],
            ['sequences.sample_id', 'sequences.ai'],
            name='seq_fkc'),
        {'mysql_row_format': 'DYNAMIC'}
    )

    sample_id = Column(Integer, autoincrement=False)
    seq_ai = Column(Integer, autoincrement=False)
    seq = relationship(Sequence, backref=backref('collapse', uselist=False))

    collapse_to_subject_sample_id = Column(Integer)
    collapse_to_subject_seq_ai = Column(Integer, index=True)
    collapse_to_subject_seq_id = Column(String(64))  # Denormalized
    instances_in_subject = Column(Integer, server_default='0', nullable=False)
    copy_number_in_subject = Column(Integer, server_default='0',
                                    nullable=False, index=True)
    samples_in_subject = Column(Integer, server_default='0', nullable=False)

    @property
    def collapse_to_seq(self):
        """Returns the sequence being collapse to"""
        return Session.object_session(self).query(Sequence).filter(
            Sequence.sample_id == self.collapse_to_subject_sample_id,
            Sequence.ai == self.collapse_to_subject_seq_ai
        ).one()


def check_string_length(cls, key, inst):
    """Checks if a string can properly fit into a given field.  If it is \
    too long, a ValueError is raised.  This prevents MySQL from truncating \
    fields that are too long.

    """
    prop = inst.prop
    # Only interested in simple columns, not relations
    if isinstance(prop, ColumnProperty) and len(prop.columns) == 1:
        col = prop.columns[0]
        # if we have string column with a length, install a length validator
        if isinstance(col.type, String) and col.type.length:
            max_length = col.type.length

            def set_(instance, value, oldvalue, initiator):
                if value is not None and len(value) > max_length:
                    msg = 'Length {} exceeds max {} for column {}'
                    raise ValueError(msg.format(len(value), max_length,
                                                col.name))
            event.listen(inst, 'set', set_)


event.listen(Base, 'attribute_instrument', check_string_length)
