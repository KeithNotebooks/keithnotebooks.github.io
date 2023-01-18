import setuptools
from setuptools import setup


with open("README.md", "r") as fh:
    long_description = fh.read()


setup(
    name="assemblypathway",
    version="0.0.1",
    author="Keith Patarroyo",
    author_email="",
    description="Interactive Plotting of Assembly Spaces",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/KeithNotebooks/keithnotebooks.github.io/",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    test_suite="test"
)
