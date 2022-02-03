from distutils.core import setup

setup(
    name='molplotly',
    packages=['molplotly'],
    version='0.1',
    # Chose a license from here: https://help.github.com/articles/licensing-a-repository
    license='Apache License, Version 2.0',

    description='TYPE YOUR DESCRIPTION HERE',
    author='William McCorkindale',
    author_email='wjm41@cam.ac.uk',

    url='https://gitlab.com/wjm41/molplotly',
    download_url='https://github.com/user/reponame/archive/v_01.tar.gz',

    keywords=['science', 'chemistry', 'cheminformatics'],
    install_requires=[
        'validators',
        'beautifulsoup4',
    ],
    classifiers=[
        
        'Development Status :: 4 - Beta',
        
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: Apache Software License'
       
        'Programming Language :: Python :: 3',
    ],
)
