from distutils.core import setup
from setuptools import setup

# read the contents of your README file
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="molplotly",
    packages=["molplotly"],
    version="1.1.1",
    license="Apache License, Version 2.0",
    description="molplotly is an add-on to plotly built on RDKit which allows 2D images of molecules to be shown in scatterplots when hovering over the datapoints.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="William McCorkindale",
    author_email="wjm41@cam.ac.uk",
    url="https://gitlab.com/wjm41/molplotly",
    download_url="https://gitlab.com/wjm41/molplotly/-/archive/v1.1.1/molplotly-v1.1.1.tar.gz",
    install_requires=["jupyter-dash", "plotly", "pandas", "ipykernel", "scikit-learn"],
    keywords=["science", "chemistry", "cheminformatics"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
    ],
)
