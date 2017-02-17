import immunedb.api.queries as queries
from immunedb.common.models import CloneStats
from immunedb.util.nested_writer import NestedCSVWriter


class MutationExporter(object):
    def __init__(self, session, clone_ids, limit_sample_ids, thresh_type,
                 thresh_value):
        self._session = session
        self._clone_ids = clone_ids
        self._limit_sample_ids = limit_sample_ids
        self._thresh_type = thresh_type
        self._thresh_value = thresh_value

        headers = ['clone_id', 'sample_id', 'total_seqs', 'region']
        self._numbers = ('unique', 'total')
        self._mtypes = ('synonymous', 'nonsynonymous', 'conservative',
                        'nonconservative', 'unknown')
        for number in self._numbers:
            for mtype in self._mtypes:
                headers.append('{}_{}'.format(mtype, number))

        self._csv = NestedCSVWriter(headers, {
            'nonsynonymous_unique': (
                lambda r: r.get('conservative_unique', 0) +
                r.get('nonconservative_unique', 0)),
            'nonsynonymous_total': (
                lambda r: r.get('conservative_total', 0) +
                r.get('nonconservative_total', 0)),
            'sample_id': lambda r: r['sample_id'] if r['sample_id'] else 'ALL'
        }, streaming=True)

    def _sample_rows(self, cid, sample_id):
        result = queries.get_clone_mutations(
            self._session, cid, self._thresh_type,
            self._thresh_value, sample_id=sample_id)

        for region in sorted(result['regions'].keys()):
            stats = result['regions'][region]
            row = {
                'clone_id': cid,
                'region': region,
                'sample_id': sample_id,
                'total_seqs': result['total_seqs'],
            }

            for number in self._numbers:
                for mtype, count in stats['counts'][number].iteritems():
                    row['{}_{}'.format(mtype, number)] = count
            self._csv.add_row(row, write_default=True, default=0,
                              write_if_stream=False)

    def get_data(self):
        for cid in sorted(self._clone_ids):
            # Write mutations for entire clone
            self._sample_rows(cid, None)
            yield self._csv.get_value()
            # Write per sample mutations
            if self._limit_sample_ids is not None:
                sample_ids = self._limit_sample_ids
            else:
                sample_ids = map(lambda r: r.sample_id, self._session.query(
                    CloneStats.sample_id).filter(
                        CloneStats.clone_id == cid,
                        CloneStats.sample_id != 0).all())
            for sample_id in sorted(sample_ids):
                self._sample_rows(cid, sample_id)
                yield self._csv.get_value()
