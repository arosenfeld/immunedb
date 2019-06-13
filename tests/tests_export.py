import io
import difflib
import os
import time
import requests
import unittest
import zipfile


def request(endpoint, query=None):
    return requests.get('http://localhost:8891' + endpoint, params=query)


class ExportTest(unittest.TestCase):
    @classmethod
    def tearDownClass(cls):
        try:
            request('/shutdown')
        except requests.exceptions.ConnectionError:
            return

    def check(self, expected_path, url, query=None):
        path = 'tests/data/export/' + expected_path
        response = request(url, query)
        assert response.status_code == 200
        uid = response.json()['uid']
        while True:
            time.sleep(1)
            resp = request('/export/job_log/' + uid)
            if resp.status_code == 200 and resp.json()['complete']:
                break

        response = request('/export/job/' + uid)
        assert response.status_code == 200

        if os.getenv('GENERATE'):
            with open(path, 'wb+') as fh:
                fh.write(response.content)
        else:
            actual_zip = zipfile.ZipFile(io.BytesIO(response.content))
            with open(path, 'rb') as fh:
                expected_zip = zipfile.ZipFile(fh)
                assert (set(expected_zip.namelist()) ==
                        set(actual_zip.namelist()))
                for name in expected_zip.namelist():
                    actual_val = actual_zip.read(name)
                    expected_val = expected_zip.read(name)
                    if expected_val != actual_val:
                        diff = difflib.unified_diff(
                            str(actual_val, encoding='utf8').split(),
                            str(expected_val, encoding='utf8').split()
                        )
                        raise AssertionError('Invalid file {}:\n{}'.format(
                            name, '\n'.join(diff)))

    def test_sequences(self):
        for schema in ('changeo', 'airr'):
            print('checking {}'.format(schema))
            self.check(
                schema + '.zip',
                '/export/sequences',
                {'format': schema}
            )

    def test_clones(self):
        for schema in ('vdjtools', 'immunedb'):
            print('checking {}'.format(schema))
            self.check(
                schema + '.zip',
                '/export/clones',
                {'format': schema}
            )

    def test_sample(self):
        self.check('sample_metadata.zip', '/export/samples')
