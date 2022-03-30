from setuptools import find_packages, setup


NAME = "aether-analysis"
DESCRIPTION='Utilities for the Analysis Pipeline'
LICENSE='Aether Biomachines Private repo'
REQUIRED = []
VERSION = 0.9

setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    packages=find_packages(),
    include_package_data=True,
    install_requires=REQUIRED,
    license=LICENSE,
    classifiers=[
        'Programming Language :: Python :: 3',
    ],
)
