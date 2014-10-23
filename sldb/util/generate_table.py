import argparse
import re
import sys
import subprocess
from tabulate import tabulate

def _get_info(path):
    date = re.search(r'\d+-\d+-\d+', path).group(0)
    ident = path.split('/')[-1].replace('.trim', '')

    return date, ident


def _make_table(paths, incl_remote):
    naming_funcs = {
        'PPG': _get_info,
        'Lupus': _get_info,
        'Vk8': None
    }
    ignores = ['AppleDouble', '3away', '3.analysis']

    rows = []
    for p in paths:
        p = p.replace('/masterTable.txt', '').strip()
        if any(e in p for e in ignores):
            continue

        for study in naming_funcs:
            if study in p:
                if naming_funcs[study] is not None:
                    date, ident = naming_funcs[study](p)
                    r = [study, ident, date]
                    if incl_remote:
                        r.append(p)
                    rows.append(r)
                break

    tbl = tabulate(rows, tablefmt='pipe').strip()
    # This monstrosity removes the leading and trailing '|' from each line
    tbl = '\n'.join(map(lambda l: l[1:-1].strip(), 
                        tbl[tbl.index('\n'):].split('\n')))

    return tbl.strip()

def run_generate_table():
    parser = argparse.ArgumentParser(
        description='Generates a study|sample|date or study|sample|date|path '
        'table from a list of paths as formatted')
    parser.add_argument('-p', action='store_true', default=False, 
                        help='Include remote path in table.')

    args = parser.parse_args()

    paths = []
    paths = sys.stdin.readlines()
    print _make_table(paths, args.p)
