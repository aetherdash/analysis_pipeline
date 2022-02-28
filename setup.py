from setuptools import find_packages, setup


NAME = "aether-analytics-utils"
DESCRIPTION='Utilities for interfacing with analytics database'
LICENSE='Aether Biomachines Private repo'
REQUIRED = []
VERSION = 0.2

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
