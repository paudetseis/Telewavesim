import setuptools
from numpy.distutils.core import setup, Extension


ext = [Extension(name='telewavesim.rmat_f',
                 sources=['src/rmat.f90', 'src/rmat_sub.f90'],
                 extra_f90_compile_args=["-O3"],
                 extra_link_args=["-L/usr/local/lib", "-lfftw3"])]

setup(
  name = 'telewavesim',
  version = '0.0.1',
  description = 'Python package for teleseismic body-wave modeling',
  author = 'Pascal Audet',
  author_email = 'pascal.audet@uottawa.ca',
  url = 'https://gitlab.com/uottawa-geophysics/SeismoPy/TeleWaveSim', 
  classifiers = [
    'Development Status :: 1 - Alpha',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Fortran',
    'Programming Language :: Python :: 3.7'],
  install_requires=["numpy>=1.16", "obspy>=1.1.0", "pyfftw>=0.11.1", "matplotlib"],
  python_requires='>=3.7',
  ext_modules = ext,
  packages=['telewavesim']
)
