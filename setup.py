from setuptools import setup, find_packages

setup(
    name="forced_phot",
    version="0.1.0",
    description="Simple forced photometry on FITS images with optional source clustering.",
    author="David Kaplan",
    author_email="kaplan@uwm.edu",
    url="",
    packages=find_packages(exclude=("tests")),
)
