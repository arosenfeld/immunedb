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
        response = request(url)
        assert response.status_code == 200
        response = response.text.strip()
        if os.getenv('GENERATE'):
            with open(path, 'w+') as fh:
                fh.write(response)
        else:
            with open(path) as fh:
                expected = fh.read().strip()
                if expected != response:
                    for e, r in zip(expected.split('\n'),
                                    response.split('\n')):
                        if e != r:
                            print('expected: ' + e)
                            print('response: ' + r)
                            break
                assert expected == response

    def test_sequences(self):
        for schema in ('changeo', 'airr'):
            print('checking {}'.format(schema))
            self.check(
                schema + '.tsv',
                '/export/sequences/' + schema
            )

    def test_clones(self):
        for schema in ('summary', 'overlap'):
            print('checking {}'.format(schema))
            self.check(
                schema + '.tsv',
                '/export/clones/' + schema
            )

    def test_sample(self):
        self.check('sample_metadata.tsv', '/export/samples')
