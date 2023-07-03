from setuptools import setup, find_packages

required = ['numpy', 'matplotlib', 'click', 'pandas>=1.1.0', 'scipy', 'biopython']

with open("README.md", "r", encoding='utf-8') as fh:
    long_description = fh.read()


setup(
    name='AKRUP',
    version='0.0.1',
    author="wangjiaqi",
    description="Ancestral Karyotype Reconstruction Universal Pipelines",
    packages=find_packages(),
    include_package_data=True,
    license="BSD License",
    long_description=long_description,
    long_description_content_type="text/markdown",
    install_requires=required,
    # package_data={'': ['*.conf', '*.pl', '*.pm','*.R', '*.exe'], 'ini': ['*']},
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
    # zip_safe=True
)
