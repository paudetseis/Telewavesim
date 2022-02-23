import os.path
import re
from numpy.distutils.core import setup, Extension
from numpy.distutils.system_info import get_info

def find_version(*paths):
    fname = os.path.join(os.path.dirname(__file__), *paths)
    with open(fname, encoding='utf-8') as fp:
        code = fp.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", code, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")


ext = [Extension(name='telewavesim.rmat_f',
                 sources=['src/rmat.f90', 'src/rmat_sub.f90'],
                 libraries=['lapack'],
                 library_dirs=get_info('lapack_opt', 1).get('library_dirs')
                 )]

setup(
    name='telewavesim',
    version=find_version('telewavesim', '__init__.py'),
    description='Python package for teleseismic body-wave modeling',
    author='Pascal Audet, Colin J. Thomson, Michael G. Bostock, Tom Eulenfeld',
    maintainer='Pascal Audet',
    author_email='pascal.audet@uottawa.ca',
    url='https://github.com/paudetseis/Telewavesim',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Fortran',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        ],
    install_requires=['numpy>=1.15', 'obspy>=1.0.0', 'matplotlib'],
    python_requires='>=3.5',
    tests_require=['pytest'],
    ext_modules=ext,
    packages=['telewavesim', 'telewavesim.tests'],
    package_data={
        'telewavesim': [
            'examples/*.ipynb',
            'examples/models/*.txt',
            'examples/Notebooks/*.ipynb']
    }
)
