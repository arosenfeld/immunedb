import sldb.api.queries as queries
from sldb.common.models import CloneStats
from sldb.common.mutations import threshold_mutations
from sldb.util.nested_writer import NestedCSVWriter


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
            'sample_id': lambda r: r['sample_id'] if r['sample_id'] else 'All'
        }, streaming=True)

    def _sample_rows(self, cid, sample_id):
        try:
            all_mutations, total_seqs = queries.get_clone_mutations(
                self._session, cid, sample_id)
        except Exception as ex:
            return

        if self._thresh_type == 'seqs':
            min_seqs = self._thresh_value
        else:
            min_seqs = int(math.ceil(self._thresh_value / 100.0 * total_seqs))

        mutations = threshold_mutations(all_mutations, min_seqs)
        for region, stats in mutations.iteritems():
            row = {
                'clone_id': cid,
                'region': region,
                'sample_id': sample_id,
                'total_seqs': total_seqs,
            }

            for number in self._numbers:
                for mtype, count in stats['counts'][number].iteritems():
                    row['{}_{}'.format(mtype, number)] = count
            self._csv.add_row(row, write_default=True, default=0,
                              write_if_stream=False)

    def get_data(self):
        for cid in self._clone_ids:
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
            for sample_id in sample_ids:
                self._sample_rows(cid, sample_id)
                yield self._csv.get_value()
