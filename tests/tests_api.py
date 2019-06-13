import json
import os
import requests
import unittest


def compare_dicts(d1, d2):
    assert set(d1.keys()) == set(d2.keys())
    for k in list(d1.keys()):
        compare_objs(d1[k], d2[k])


def compare_objs(o1, o2):
    assert (type(o1) == type(o2) or
            (type(o1) == type(o2) == str))
    if type(o1) == dict:
        compare_dicts(o1, o2)
    elif type(o1) == list:
        assert len(o1) == len(o2)
        for item1, item2 in zip(o1, o2):
            compare_objs(item1, item2)
    else:
        assert o1 == o2, (o1, o2)


def pp(j):
    return json.dumps(j, sort_keys=True, indent=4, separators=(',', ': '))


class ApiTest(unittest.TestCase):
    def request(self, endpoint, data={}):
        return requests.post('http://localhost:8891' + endpoint, json=data)

    def check(self, expected_path, url, data={}):
        path = 'tests/data/responses/' + expected_path + '.json'
        response = self.request(url, data).json()
        if os.getenv('GENERATE'):
            with open(path, 'w+') as fh:
                fh.write(pp(response))
        else:
            with open(path) as fh:
                expected = json.loads(fh.read())
            try:
                compare_objs(expected, response)
            except AssertionError:
                with open('assert.log', 'w+') as err:
                    err.write('Expected\n{}\n\nResponse\n{}'.format(
                        pp(expected), pp(response)))
                raise

    def test_endpoints(self):
        endpoints = {
            'sample_list': '/samples/list',
            'sequences_list': '/sequences/list',
            'sequence_check':
                '/sequence/2/M03592:18:000000000-AHJWK:1:1102:9856:7192',

            'clones_list': '/clones/list',
            'clones_list_subject1': '/clones/list/1',
            'clone': '/clone/10',

            'analyze': '/samples/analyze/T2',
            'overlap': '/samples/overlap/T2',
            'v_usage': '/samples/v_usage/T2',

            'subject_list': '/subjects/list',
            'subject': '/subject/1',
        }
        if os.environ.get('BASELINE_PATH'):
            endpoints['clone_pressure'] = '/clone/pressure/1'
        for cid in range(1, 25):
            for value in ('sequences', 'mutations'):
                name = ('clone', value, str(cid))
                endpoints['_'.join(name)] = '/' + '/'.join(name)

        for check, endpoint in endpoints.items():
            print(check, endpoint)
            self.check(check, endpoint)

        self.check('clones_filters1', '/clones/list', {
            'filters': {
                'min_cdr3_num_nts': 5,
                'max_cdr3_num_nts': 30,
                'min_size': 2,
                'max_size': 100,
                'size_field': 'uniques'
            }
        })
        self.check('clones_filters2', '/clones/list', {
            'filters': {
                'subject_id': 1,
                'id': 1
            }
        })
        self.check('clones_filters3', '/clones/list', {
            'order_field': 'id',
            'order_dir': 'asc'
        })
        self.check('clones_filters4', '/clones/list', {
            'order_field': 'id',
            'order_dir': 'desc'
        })
