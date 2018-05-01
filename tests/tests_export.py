import os
import requests
import unittest


def request(endpoint):
    return requests.get('http://localhost:8891' + endpoint)


class ExportTest(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        request('/shutdown')

    def check(self, expected_path, url):
        path = 'tests/data/export/' + expected_path
        response = request(url).text
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
        for schema in ('changeo', 'airr'):
            print 'checking {}'.format(schema)
            self.check(
                schema + '.tsv',
                '/export/tsv/' + schema
            )
