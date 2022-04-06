from setuptools import find_packages, setup

setup(
    name="molplotly",
    version="1.1.2",
    description="plotly add-on to render molecule images on mouseover",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/wjm41/molplotly",
    author="William McCorkindale",
    license="Apache License 2.0",
    packages=find_packages(),
    install_requires=[
        "dash==2.0.0",
        "werkzeug==2.0.0",
        "pytest-cov ~= 3.0.0",
        "jupyter-dash",
        "plotly~=5.6.0",
        "pandas~=1.4.1",
        "ipykernel~=6.9.1",
        "scikit-learn~=1.0.0",
    ],
    extras_require={"test": ["pytest"]},
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering",
    ],
)
