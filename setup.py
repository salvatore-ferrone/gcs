from setuptools import setup, find_packages

setup(
    name='gcs',
    version='0.1',
    packages=find_packages(),
    description='A simple package for data extraction and path handling in GCS projects',
    author='Salvatore Ferrone',
    author_email='salvatore.ferrone@uniroma1.it',
    keywords='gcs data extraction path handling',
    install_requires=[
        "numpy"
    ],
    package_data={
        # If any package contains *.yaml files, include them:
        'gcs': ['*.yaml'],
        # And if you need to be more specific, you can specify the package name:
        # 'gcs': ['paths.yaml'],
    },
)