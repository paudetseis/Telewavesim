import setuptools
from numpy.distutils.core import setup, Extension

version=open('version.txt').read().split()[0]

ext = [Extension(name='telewavesim.rmat_f',
                sources=['src/rmat.f90', 'src/rmat_sub.f90'],
                extra_f90_compile_args=["-O3"],
                # extra_link_args=["-lfftw3"])]
                # For separate Fortran build, fftw and lapack libraries,
                # change the location of the Library path where
                # the `fftw3` and `lapack` libraries reside and uncomment the
                # following line. Don't forget to 
                # comment out the previous line above.
                extra_link_args=["-L/usr/local/lib", "-lfftw3", "-llapack"])]

setup(
  name = 'telewavesim',
  version = version,
  description = 'Python package for teleseismic body-wave modeling',
  author = 'Pascal Audet, Colin J. Thomson, Michael G. Bostock',
  maintainer = 'Pascal Audet',
  author_email = 'pascal.audet@uottawa.ca',
  url = 'https://gitlab.com/uottawa-geophysics/SeismoPy/TeleWaveSim', 
  classifiers = [
    'Development Status :: 0 - Alpha',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Fortran',
    'Programming Language :: Python :: 3.7'],
  install_requires = ["numpy>=1.15", "obspy>=1.1.0", "pyfftw>=0.11.1", "matplotlib"],
  python_requires = '>=3.7',
  ext_modules = ext,
  packages = ['telewavesim'],
  package_data = {
    "telewavesim": [
      "examples/*.ipynb",
      "examples/models/*.txt",
      "examples/Notebooks/*.ipynb",
      "examples/Notebooks/*.ipynb"]
  }
)
