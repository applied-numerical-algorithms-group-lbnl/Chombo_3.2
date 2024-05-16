
from distutils.core import setup

setup(
    name='AMRFile',
    version='0.1.0',
    author='Stephen l Cornford',
    author_email='ggslc@bristol.ac.uk',
    packages=['amrfile'],
    url='http://bisicles.lbl.gov',
    license='LICENSE.txt',
    description='Access to Chombo AMR data stored in hdf5 files',
    long_description=open('README.txt').read(),
    install_requires=[
        "numpy >= 1.1.0",
    ],
)
