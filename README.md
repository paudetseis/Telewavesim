# Telewavesim: Software for teleseismic body wave modeling

![](examples/picture/tws_logo.png)

The structure of the Earth's crust and upper mantle gives useful information on the 
internal composition and dynamics of our planet. Some of the most widely used techniques
to infer these properties are based on examining the effect of teleseismic body wave 
(i.e., P and S waves that originate from distant earthquakes and arrive as plane waves)
propagation (e.g., transmission and scattering) through stratified media. Modeling the 
seismic response from stacks of subsurface layers is therefore an essential tool in 
characterizing their effect on observed seismograms.

This package contains `python` and `fortran` modules to synthesize teleseismic 
body-wave propagation through stacks of generally anisotropic and strictly horizontal 
layers using the matrix propagator approach of [Kennett (1983)](#references). 
The software also properly models reverberations from an overlying column of water, 
effectively simulating ocean-bottom seismic (OBS) station recordings. The software 
will be useful in a variety of teleseismic receiver-based studies, such as P or S 
receiver functions, long-period P-wave polarization, shear-wave splitting from 
core-refracted shear waves (i.e., SKS, SKKS), etc. It may also be the starting point 
for stochastic inverse methods (e.g., Monte Carlo sampling). The main part of the
code is written in `fortran` with `python` wrappers. Common computational 
workflows are covered in the Jupyter notebooks bundled with this package.

## Navigation / Notebooks

Included in this package is a set of Jupyter Notebooks, which give examples on how to call the various routines and obtain plane wave seismograms and receiver functions. The Notebooks desribe how to reproduce published examples of synthetic data from [Audet (2016)](#references) and [Porter et al. (2011)](#references).

- [sim_obs_Audet2016.ipynb](./examples/Notebooks/sim_obs_Audet2016.ipynb): Example plane wave seismograms and P receiver functions for OBS data from [Audet (2016)](#Audet).
- [sim_Prfs_Porter2011.ipynb](./examples/Notebooks/sim_Prfs_Porter2011.ipynb): Example P receiver functions from [Porter et al. (2011)](#Porter)
- [sim_SKS.ipynb](./examples/Notebooks/sim_SKS.ipynb): Example plane wave seismograms for SKS splitting studies.

## Dependencies

The current version was developed using **Python3.7**
Also, the following packages are required:

- [`obspy`](https://github.com/obspy/obspy/wiki)
- [`pyfftw`](https://pyfftw.readthedocs.io/en/latest/)
- [`fftw`](http://www.fftw.org)
- [`lapack`](http://www.netlib.org/lapack)
- A working `Fortran` compiler (e.g., GCC (gfortran))

By  default, both `numpy` and `matplotlib` are installed as dependencies of `obspy`. 
See below for full installation details.

## Installation

There is more than one way to install the software. We highly recommend installing 
the software and its dependencies using the `conda` package manager in a virtual environment. 
However, we recognize that some users may prefer to use a different `Fortran` compiler 
(e.g., Intel (c) Fortran), and we provide steps below to install the software using 
a system-wide available `Fortran` compiler.

Regardless of your choice, start with this step:

- Clone the repository:
```bash
$ git clone https://github.com/paudetseis/Telewavesim.git
$ cd Telewavesim
```

### Conda installation

- Add the `conda-forge` channel:
```bash
$ conda config --add channels conda-forge
```

- We advise creating a custom 
[conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
where `telewavesim` can be installed along with its dependencies. 

* Create a new environment and install all dependencies:
```bash
$ conda create -n tws python=3.7 obspy pyfftw
```
* or create it from the `tws_env.yml` file:
```bash
$ conda env create -f tws_env.yml
```
* Activate the newly created environment:
```bash
$ conda activate tws
(tws) $
```

Here, depending on your preference, you can further use `conda` to install the required Fortran 
compiler and the `fftw` library (the `lapack` library will be already installed from the 
previous steps) or use a pre-existing `Fortran` compiler 
(see [Separate Fortran build](#separatefortranbuild)). Using a `conda` environment will
will point to the correct path for dynamic linking. On a MacOSX, the `gfortran` package is `gfortran_osx-64`; for Linux, the `gfortran` package is `gfortran_linux-64` (check out https://anaconda.org/search?q=gfortran for the available packages).

- Install `gfortran` and the `fftw` library:

```bash
(tws) $ conda install gfortran_osx-64 fftw
```

You can check that the active Fortran compiler is the `conda` version residing in 
the `tws` path:
```bash
(tws) $ which gfortran
```

Install the software using `pip`, or build and install from source

- Using `pip`:
```bash
(tws) $ pip install .
```
- Build and install the project:

```bash
(tws) $ python setup.py build 
(tws) $ python setup.py install
```

### Separate Fortran build

For `conda` users, we recommend creating and activating the `conda` environment as 
detailed above. 

- You will need to use a system-wide available Fortran compiler (e.g., located in 
`/usr/local/bin`), and independently download and install the `fftw` library 
(http://www.fftw.org).

- Once this is done, edit the `setup.py` file to modify the content of `extra_link_args` in the `Extension` class to point to the location of your compiled `fftw` and `lapack` libraries. In the example below the libraries are installed in `/usr/local/lib`:

```python
ext = [Extension(name='telewavesim.rmat_f',
                 sources=['src/rmat.f90', 'src/rmat_sub.f90'],
                 extra_f90_compile_args=["-O3"],
                 extra_link_args=["-L/usr/local/lib", "-lfftw3", "-llapack"])]
```

Install the software using `pip`, or build and install from source

- Using `pip`:
```bash
$ pip install .
```
- Build and install the project:

```bash
$ python setup.py build 
$ python setup.py install
```

---
**NOTE**

If you are actively working on the code, or making frequent edits, it is advisable to perform 
the ``pip`` installation with the `-e` flag: 
```bash
pip install -e .
```

This enables an editable installation, where symbolic links are used rather than straight 
copies. This means that any changes made in the local folders will be reflected in the 
package available on the system.

---

## References
1. Audet, P. (2016). Receiver functions using OBS data: promises and limitations from numerical modelling and examples from the Cascadia Initiative. Geophysical Journal International, 205, 1740-1755. https://doi.org/10.1093/gji/ggw111
2. Kennett, B.L.N. (1983). Seismic wave propagation in stratified media. Cambridge University Press, 342pp.
3. Porter, R., Zandt, G., & McQuarrie, N. (2011). Pervasive lower-crustal seismic anisotropy in Southern California: Evidence for underplated schists and active tectonics. Lithosphere, 3(3), 201-220. https://doi.org/10.1130/L126.1
