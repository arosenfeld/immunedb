import os

from setuptools import setup, Extension

if os.environ.get('READTHEDOCS'):
    install_requires = ['sphinxcontrib-websupport']
    ext_modules = []
else:
    with open('requirements.txt') as req:
        install_requires = [l.strip() for l in req]
    ext_modules = [
        Extension('dnautils', sources=['lib/dnautils.c'],
                  extra_compile_args=['-std=c99'])
    ]

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()

setup(
    name='ImmuneDB',
    version='0.28.2',
    author='Aaron M. Rosenfeld',
    author_email='ar374@drexel.edu',
    url='https://github.com/arosenfeld/immunedb',
    packages=[
        'immunedb',
        'immunedb.aggregation',
        'immunedb.api',
        'immunedb.common',
        'immunedb.exporting',
        'immunedb.exporting.clones',
        'immunedb.identification',
        'immunedb.importing',
        'immunedb.trees',
        'immunedb.util',
    ],
    scripts=[
        'bin/immunedb_admin',
        'bin/immunedb_clones',
        'bin/immunedb_clone_import',
        'bin/immunedb_clone_stats',
        'bin/immunedb_clone_pressure',
        'bin/immunedb_clone_trees',
        'bin/immunedb_collapse',
        'bin/immunedb_export',
        'bin/immunedb_local_align',
        'bin/immunedb_identify',
        'bin/immunedb_import',
        'bin/immunedb_metadata',
        'bin/immunedb_modify',
        'bin/immunedb_rest',
        'bin/immunedb_sample_stats',
        'bin/immunedb_sql',
        'bin/immunedb_genotype',
        'bin/run_tigger'
    ],
    install_requires=install_requires,
    ext_modules=ext_modules,
    license='LICENSE.txt',
    description='''A system for the analysis and exploration of high-throughput
        adaptive immune receptor sequencing data''',
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
