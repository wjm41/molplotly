from distutils.core import setup

setup(
    name='molplotly',
    packages=['molplotly'],
    version='0.11',
    license='Apache License, Version 2.0',

    description='molplotly is an add-on to plotly built on RDKit which allows 2D images of molecules to be shown in scatterplots when hovering over the datapoints.',
    author='William McCorkindale',
    author_email='wjm41@cam.ac.uk',

    url='https://gitlab.com/wjm41/molplotly',
    download_url='https://gitlab.com/wjm41/molplotly/-/archive/v_011/molplotly-v_011.tar.gz',

    keywords=['science', 'chemistry', 'cheminformatics'],
    install_requires=[
        'jupyter-dash',
        'pandas'
    ],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
    ]
)
