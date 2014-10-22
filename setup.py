from setuptools import setup

setup(name='SLDB',
    version='1.0.0',
    author='Aaron M. Rosenfeld',
    author_email='ar374@drexel.edu',
    packages=[
            'sldb',
    ],
    scripts=['bin/mt2db'],
    install_requires=[
        'sqlalchemy>=0.9.8',
        'flask',
        'flask-restless',
    ],
    license='LICENSE.txt',
    description='Various utilities for Drexel\'s Systems Immunology Lab.',
)
