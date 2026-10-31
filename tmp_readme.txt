# just for running stuff on legendre:

# e.g. to compute generators for degree 3:
PYTHONPATH=~/sage-extra sage -python scripts/find_generators.py --degrees 3 --max-disc 10^6 --only-new

# to run tests:
sage -python -m unittest discover tests -v
