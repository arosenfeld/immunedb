from setuptools import setup

setup(name='SLDB',
      version='1.0.5',
      author='Aaron M. Rosenfeld',
      author_email='ar374@drexel.edu',
      packages=[
              'sldb',
              'sldb.api',
              'sldb.common',
              'sldb.conversion',
              'sldb.trees',
              'sldb.util',
      ],
      scripts=[
          'bin/sldb_mt2db',
          'bin/sldb_rest',
          'bin/sldb_newick2json'
      ],
      install_requires=[
          'sqlalchemy>=0.9.8',
          'Flask',
          'Flask-SQLAlchemy',
          'flask-restless',
          'ete2 >= 2.2',
      ],
      license='LICENSE.txt',
      description='Various utilities for Drexel\'s Systems Immunology Lab.')
