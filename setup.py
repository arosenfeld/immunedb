from setuptools import setup, Extension

dnautils = Extension('dnautils', sources=['lib/dnautils.c'],
                     extra_compile_args=['-std=c99'])

setup(
    name='SLDB',
    version='0.12.3',
    author='Aaron M. Rosenfeld',
    author_email='ar374@drexel.edu',
    packages=[
        'sldb',
        'sldb.aggregation',
        'sldb.api',
        'sldb.common',
        'sldb.exporting',
        'sldb.identification',
        'sldb.importing',
        'sldb.trees',
        'sldb.util',
    ],
    scripts=[
        'bin/sldb_admin',
        'bin/sldb_clones',
        'bin/sldb_clone_stats',
        'bin/sldb_clone_pressure',
        'bin/sldb_clone_trees',
        'bin/sldb_collapse',
        'bin/sldb_local_align',
        'bin/sldb_identify',
        'bin/sldb_import_delim',
        'bin/sldb_hvquest',
        'bin/sldb_rest',
        'bin/sldb_sample_stats',
    ],
    install_requires=[
        'gevent',
        'sqlalchemy',
        'biopython',
        'bottle',
        'ete2',
        'numpy',
        'scipy',
        'PyMySQL',
    ],
    ext_modules=[dnautils],
    license='LICENSE.txt',
    description='Various utilities for Drexel\'s Systems Immunology Lab.'
)
