## Telewavesim: Python software for teleseismic body wave modeling

This python package contains modules to synthesize teleseismic (i.e., plane-wave) 
body-wave propagation through stacks of generally anisotropic layers using the
matrix propagator approach of Kennett (1983). The main subroutine Rmatrix
was originally written by Colin Thomson (c). The software also properly models water 
column reverberations, simulating ocean-bottom seismic (OBS) station recordings. The software
will be useful in a variety of receiver-based studies, such as P or S receiver functions,
long-period P-wave polarization, shear-wave splitting from core-refracted shear waves (i.e., SKS, SKKS),
etc. 


Pascal Audet, August 2019



## Installation

### Dependencies

You will need **Python 3.7+**.
Also, the following packages are required:

- [`numpy`](http://numpy.org)
- [`matplotlib`](https://matplotlib.org/)
- [`obspy`](https://github.com/obspy/obspy/wiki)
- [`pyfftw`](https://pyfftw.readthedocs.io/en/latest/)
- [`fftw3`](http://www.fftw.org)

Note that both `numpy` and `matplotlib` are installed as dependencies of `obspy`. See below for full installation details. You also need to download and install the [`fftw3`](http://www.fftw.org) library independently. 

### Conda environment

You can create a custom [conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
where `telewavesim` can be installed along with its dependencies.

Clone the repository:
```bash
git clone https://gitlab.com/uottawa-geophysics/SeismoPy/TeleWaveSim.git
cd telewavesim
```

Create the environment from the `tws_env.yml` file:
```bash
conda env create -f tws_env.yml
```
Activate the newly created environment:
```bash
conda activate tws
```

Finally, edit the setup.py file to modify the link arguments in the Extension class to point to your compiled `fftw3` library:

```python
ext = [Extension(name='telewavesim.rmat_f',
                 sources=['src/rmat.f90', 'src/rmat_sub.f90'],
                 extra_f90_compile_args=["-O3"],
                 extra_link_args=["-L/usr/local/lib", "-lfftw3"])]
```

### 1) Installing using pip

Once the previous steps are performed, you can install `telewavesim` using pip:
```bash
pip install .
```

### 2) Building and Installing

Alternatively, you can build and install the project (from the root of the source tree, e.g., inside the cloned `telewavesim` directory):

```bash
python setup.py build 
python setup.py install
```

Please note, if you are actively working on the code, or making frequent edits, it is advisable
to perform the pip installation with the -e flag. This enables an editable installation, where
symbolic links are used rather than straight copies. This means that any changes made in the
local folders will be reflected in the packages available on the system.

<!-- ### Installing using pip

You can install `telewavesim` using the
[`pip package manager`](https://pypi.org/project/pip/):

```bash
pip install telewavesim
```
All the dependencies will be automatically installed by `pip`.

### Installing with conda

You can install `telewavesim` using the [conda package manager](https://conda.io).
Its required dependencies can be easily installed with:

```bash
conda install -c conda-forge obspy
conda install numpy=1.16
conda install pyfftw
```

Then `telewavesim` can be installed with `pip`:

```bash
pip install telewavesim
```

#### Conda environment

Alternatively, you can create a custom
[conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
where `telewavesim` can be installed along with its dependencies.

Clone the repository:
```bash
git clone https://github.com/brmather/pycurious
cd telewavesim
```

Create the environment:
```bash
conda create --name tws
```

Activate the newly created environment:
```bash
conda activate tws
```

And install `telewavesim` with `pip`:
```bash
pip install telewavesim
```
 -->

<!-- Both implementations:

* $ conda --add channels conda-forge
* $ conda create --name tws
* $ conda activate tws
* $ conda install obspy
* $ conda install pyfftw
* $ conda install numpy=1.16

With FORTRAN wrapper:

* Install FFTW3 library: http://www.fftw.org
* Install LAPACK library: http://www.netlib.org/lapack/
* Modify the makefile in rmat_f to point directories to libraries and include files
    * (currently in /usr/local/lib and /usr/local/include)
* $ cd rmat_f
* $ make


Package Contents:
-------------------

telewavesim Module:

* cfg:
    * conf.py -> defines global variables
* elast:
    * elast_stifness.py -> defines stifness matrices for various minerals and rocks
    * elast_util.py -> loads stifness matrices to interface with code
* green:
    * green.py -> contains functions to generate Green's functions using rmatrix for land or ocean-bottom stations
* plotting:
    * wiggle.py -> contains functions to plot synthetic traces
* rmat:
    * rmatrix.py -> matrix propagator algorithm at single frequency
    * rtfluid.py -> R/T coefficients for fluid-solid interface
* rmat_f:
    * makefile
    * rmat.f90 -> contains conf, rmat and green modules in FORTRAN
    * lib:
        * makefile
        * rmat_sub.f90 -> contains subroutines used in rmat
        * disp.f90 -> subroutines to display matrices
* utils:
    * model_util.py -> utility module to read and define model parameters
    * trace_util.py -> utility module to handle traces

examples: Python scripts making use of the module for manipulation and creation of databases

 -->