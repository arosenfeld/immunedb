import os
import requests
import unittest

from immunedb.exporting.clone_export import CloneExport
from immunedb.exporting.sequence_export import SequenceExport


def request(endpoint, data={}):
    return requests.post('http://localhost:8891' + endpoint, data=data)


class ExportTest(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        request('/shutdown')

    def check(self, expected_path, url, data={}):
        path = 'tests/data/export/' + expected_path
        response = request(url, data).text
        if os.getenv('GENERATE'):
            with open(path, 'w+') as fh:
                fh.write(response)
        else:
            with open(path) as fh:
                expected = fh.read().strip()
                if expected != response.strip():
                    print 'EXPECTED'
                    print expected
                    print 'RESPONSE'
                    print response.strip()
                assert expected == response.strip()

    def test_sequences(self):
        for ext in ('fasta', 'fastq', 'clip'):
            print 'checking {}'.format(ext)
            self.check(
                'sequences.' + ext,
                '/export/sequences/sample/T2', data={
                    'format': ext,
                    'fields': ','.join(sorted(
                        SequenceExport.allowed_fields.keys()))
                })

    def test_clones(self):
        self.check(
            'clones.csv',
            '/export/clones/sample/T2', data={
                'fields':
                ','.join(sorted(
                    set(CloneExport.allowed_fields.keys()) -
                    set(['tree'])))
            })

    def test_mutations(self):
        self.check(
            'mutations.csv',
            '/export/mutations/sample/T2', data={
                'fields': ','.join(sorted(CloneExport.allowed_fields.keys())),
            })
