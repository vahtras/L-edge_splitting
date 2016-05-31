from setuptools import setup

setup(name="Ledger",
    version="0.1",
    packages=["ledges"],
    scripts=[
        'scripts/ledges',
        'scripts/generate_spectrum.sh',
        'scripts/justplot2.py'
        ],
    author="Olav Vahtras",
    author_email="vahtras@kth.se"
    )
