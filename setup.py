from setuptools import setup, Extension

dnautils = Extension('dnautils', sources=['lib/dnautils.c'],
                     extra_compile_args=['-std=c99'])

setup(
    name='airrdb',
    version='0.17.0',
    author='Aaron M. Rosenfeld',
    author_email='ar374@drexel.edu',
    url='https://github.com/arosenfeld/airrdb',
    packages=[
        'airrdb',
        'airrdb.aggregation',
        'airrdb.api',
        'airrdb.common',
        'airrdb.exporting',
        'airrdb.identification',
        'airrdb.importing',
        'airrdb.trees',
        'airrdb.util',
    ],
    scripts=[
        'bin/airrdb_admin',
        'bin/airrdb_clones',
        'bin/airrdb_clone_stats',
        'bin/airrdb_clone_pressure',
        'bin/airrdb_clone_trees',
        'bin/airrdb_collapse',
        'bin/airrdb_local_align',
        'bin/airrdb_identify',
        'bin/airrdb_import',
        'bin/airrdb_rest',
        'bin/airrdb_sample_stats',
        'bin/airrdb_sql',
    ],
    install_requires=[
        'gevent',
        'sqlalchemy',
        'biopython',
        'bottle',
        'ete2',
        'numpy',
        'PyMySQL',
    ],
    ext_modules=[dnautils],
    license='LICENSE.txt',
    description='A module for efficient storage and analysis of '
                'high-throughput B-cell sequence data.',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'Natural Language :: English',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Bio-Informatics'
    ]
)
