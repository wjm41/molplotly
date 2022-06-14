from setuptools import find_packages, setup

setup(
    name="molplotly",
    version="1.1.4",
    description="plotly add-on to render molecule images on mouseover",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wjm41/molplotly",
    author="William McCorkindale",
    license="Apache License 2.0",
    packages=find_packages(),
    install_requires=[
        "dash>=2.0.0",
        "werkzeug>=2.0.0",
        "jupyter-dash>=0.4.2",
        "plotly>=5.0.0",
        "pandas",
        "ipykernel",
        "nbformat",
    ],
    extras_require={"test": ["pytest", "pytest-cov"]},
    keywords=["science", "chemistry", "cheminformatics"],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: Apache Software License",
        "Programming Language :: Python :: 3",
    ],
)
