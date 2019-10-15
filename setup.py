import os.path
import platform
import re
from numpy.distutils.core import setup, Extension


def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname) as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


if platform.system() == 'Linux':
    # this hack hopefully  will become unecessary in the future
    # see issue #2
    # ci fail: https://travis-ci.org/trichter/Telewavesim/builds/594469551
    from numpy import __config__
    attrs = ['blas_mkl_info', 'blas_opt_info', 'lapack_mkl_info',
             'lapack_opt_info', 'openblas_lapack_info']
    libdirs = set()
    for attr in attrs:
        obj = getattr(__config__, attr, {}).get('library_dirs')
        if obj is not None:
            libdirs.update(obj)
    ext_kw = dict(library_dirs=list(libdirs), libraries=['lapack'])
else:
    ext_kw = {}

ext = [Extension(name='telewavesim.rmat_f',
                 sources=['src/rmat.f90', 'src/rmat_sub.f90'],
                 **ext_kw)]

setup(
    name                = 'telewavesim',
    version             = find_version('telewavesim', '__init__.py'),
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
    packages            = ['telewavesim', 'telewavesim.tests'],
    package_data        = {
        'telewavesim': [
            'examples/*.ipynb',
            'examples/models/*.txt',
            'examples/Notebooks/*.ipynb']
    }
)
