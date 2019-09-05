import setuptools
from numpy.distutils.core import setup, Extension

ext = [Extension(name='telewavesim.rmat_f',
                sources=['src/rmat.f90', 'src/rmat_sub.f90'])]

setup(
    name                = 'telewavesim',
    version             = '0.1.0',
    description         = 'Python package for teleseismic body-wave modeling',
    author              = 'Pascal Audet, Colin J. Thomson, Michael G. Bostock',
    maintainer          = 'Pascal Audet',
    author_email        = 'pascal.audet@uottawa.ca',
    url                 = 'https://github.com/paudetseis/Telewavesim', 
    download_url        = 'https://github.com/paudetseis/Telewavesim/archive/0.1.0.tar.gz',
    classifiers         = [
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'],
    install_requires    = ['numpy>=1.15', 'obspy>=1.0.0', 'pyfftw>=0.11.1', 'matplotlib'],
    python_requires     = '>=3.5',
    tests_require       = ['pytest'],
    ext_modules         = ext,
    packages            = ['telewavesim'],
    package_data        = {
        'telewavesim': [
            'examples/*.ipynb',
            'examples/models/*.txt',
            'examples/Notebooks/*.ipynb']
    }
)
