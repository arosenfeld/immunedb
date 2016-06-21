import multiprocessing as mp
import unittest

from immunedb.common.config import get_base_arg_parser


class ParserTest(unittest.TestCase):
    def test_parser(self):
        desc = 'this is the description'
        parser = get_base_arg_parser(desc)
        assert parser.description == desc
        args = parser.parse_args(['config/path.json'])
        assert set(vars(args).keys()) == set(('db_config', 'nproc'))
        assert args.db_config == 'config/path.json'
        try:
            num_cpu = mp.cpu_count()
        except NotImplementedError:
            num_cpu = 4
        assert args.nproc == num_cpu
