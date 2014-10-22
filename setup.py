from setuptools import setup

setup(name='SLDB',
    version='1.0.0',
    author='Aaron M. Rosenfeld',
    author_email='ar374@drexel.edu',
    packages=[
            'sldb',
            'sldb.trees',
    ],
    scripts=[
        'bin/sldb_mt2db',
        'bin/sldb_api',
        'bin/sldb_newick2json'
    ],
    install_requires=[
        'sqlalchemy>=0.9.8',
        'Flask',
        'flask-restless',
        'ete2 >= 2.2',
    ],
    license='LICENSE.txt',
    description='Various utilities for Drexel\'s Systems Immunology Lab.',
)
