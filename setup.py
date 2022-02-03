from distutils.core import setup

setup(
    name='molplotly',         # How you named your package folder (MyLib)
    packages=['molplotly'],   # Chose the same as "name"
    version='0.1',      # Start with a small number and increase it with every change you make
    # Chose a license from here: https://help.github.com/articles/licensing-a-repository
    license='MIT',
    # Give a short description about your library
    description='TYPE YOUR DESCRIPTION HERE',
    author='William McCorkindale',                   # Type in your name
    author_email='wjm41@cam.ac.uk',      # Type in your E-Mail
    # Provide either the link to your github or to your website
    url='https://github.com/user/wjm41',
    # I explain this later on
    download_url='https://github.com/user/reponame/archive/v_01.tar.gz',
    # Keywords that define your package best
    keywords=['science', 'chemistry', 'cheminformatics'],
    install_requires=[            # I get to this in a second
        'validators',
        'beautifulsoup4',
    ],
    classifiers=[
        # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Development Status :: 4 - Beta',
        # Define that your audience are developers
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',   # Again, pick a license
        # Specify which pyhton versions that you want to support
        'Programming Language :: Python :: 3',
    ],
)
