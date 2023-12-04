import distutils.core
from distutils.core import setup
import setuptools
from setuptools import find_packages
import os

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
     name='ElastiCouplings',
     version='1.0',
     packages=['ElastiCouplings'] ,
     author="Dario Fiore Mosca",
     author_email="dario.fiore.mosca@univie.ac.at",
     description="   ",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url=" ",
     install_requires=['numpy','pandas','scipy'],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
 )
