import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="rnaseq_barcode",
    version="0.0.1",
    author="Tom Roeschinger, Manuel Razo, Scott Saunders, Rob Phillips",
    author_email="troeschi {at} caltech {dot} edu",
    description="This repository contains all active research materials for the project testing the use of RNA-seq as a tool for biophysics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mrazomej/rnaseq_barcode.git",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
