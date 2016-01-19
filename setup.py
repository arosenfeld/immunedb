from setuptools import setup, Extension

dnautils = Extension('dnautils', sources=['lib/dnautils.c'],
                     extra_compile_args=['-std=c99'])

setup(
    name='SLDB',
    version='0.15.0',
    author='Aaron M. Rosenfeld',
    author_email='ar374@drexel.edu',
    url='https://github.com/arosenfeld/sldb',
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
        'bin/sldb_import',
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
