from setuptools import setup, find_packages

required = ['numpy', 'matplotlib', 'click', 'pandas>=1.1.0', 'scipy', 'biopython']

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()


setup(
    name='AKRUP',
    version='1.0.4',
    author="wangjiaqi",
    description="Ancestral Karyotype Reconstruction Universal Pipeline",
    url="https://github.com/Genome-structure-evolution-analysis/AKRUP",
    packages=find_packages(),
    include_package_data=True,
    license="BSD License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=required,
    package_data={'': ['*.conf', '*']},
    entry_points={
        'console_scripts': [
        'AKRUP = AKRUP.run:main'
        ]
        },
    classifiers=[
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: BSD License",
        "Operating System :: OS Independent",
    ]
)
