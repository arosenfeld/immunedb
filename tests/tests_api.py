import json
import os
import requests
import unittest

import immunedb.common.config as config


def compare_dicts(d1, d2):
    assert set(d1.keys()) == set(d2.keys())
    for k in d1.keys():
        compare_objs(d1[k], d2[k])


def compare_objs(o1, o2):
    assert (type(o1) == type(o2) or
            (type(o1) in (str, unicode) and type(o2) in (str, unicode)))
    if type(o1) == dict:
        compare_dicts(o1, o2)
    elif type(o1) == list:
        assert len(o1) == len(o2)
        for item1, item2 in zip(o1, o2):
            compare_objs(item1, item2)
    else:
        assert o1 == o2


class ApiTest(unittest.TestCase):
    def tearDown(self):
        self.request('/shutdown')

    def request(self, endpoint, data={}):
        return requests.post('http://localhost:8891' + endpoint, data)

    def check(self, expected_path, url, data={}):
        path = 'tests/data/responses/' + expected_path + '.json'
        response = self.request(url, {}).json()
        if os.getenv('GENERATE'):
            with open(path, 'w+') as fh:
                json.dump(response, fh, sort_keys=True,
                          indent=4, separators=(',', ': '))
        else:
            with open(path) as fh:
                compare_objs(json.loads(fh.read()), response)

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
        for cid in range(1, 25):
            # TODO: pressure and lineages
            for value in ('sequences', 'mutations'):
                name = ('clone', value, str(cid))
                endpoints['_'.join(name)] = '/' + '/'.join(name)
        for check, endpoint in endpoints.iteritems():
            self.check(check, endpoint)
