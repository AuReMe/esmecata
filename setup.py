from setuptools import setup

setup(
    name='EsMeCaTa',
    url='https://github.com/ArnaudBelcour/esmecata',
    license='GPLv3+',
    description=
    'EsMeCaTa: Estimating Metabolic Capabilties from Taxonomy',
    author='AuReMe',
    author_email='gem-aureme@inria.fr',
    packages=['esmecata'],
    package_dir={'esmecata': 'esmecata'},
    entry_points={
        'console_scripts': [
            'esmecata = esmecata.__main__:main',
        ]
    },
    install_requires=['biopython', 'pandas', 'requests', 'ete3'],
)