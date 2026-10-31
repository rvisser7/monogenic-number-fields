# just for running stuff on legendre:

# e.g. to compute generators for degree 3:
PYTHONPATH=~/sage-extra sage -python scripts/find_generators.py --degrees 3 --max-disc 10^6 --only-new

# to run tests:
sage -python -m unittest discover tests -v

# to check/verify the data:
sage -python scripts/check_data.py

# to check/verify data using lmfdb:
PYTHONPATH=~/sage-extra sage -python scripts/check_data.py --source db
