# Telewavesim: Software for teleseismic body wave modeling

![](./telewavesim/examples/picture/tws_logo.png)

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

[![DOI](https://zenodo.org/badge/204565459.svg)](https://zenodo.org/badge/latestdoi/204565459)
[![PyPI version](https://badge.fury.io/py/telewavesim.svg)](https://badge.fury.io/py/telewavesim)

## Installation

### Dependencies

The current version was developed using **Python3.7**
Also, the following packages are required:

- [`gfortran`](https://gcc.gnu.org/wiki/GFortran) (or any Fortran compiler)
- [`obspy`](https://github.com/obspy/obspy/wiki)
- [`pyfftw`](https://pyfftw.readthedocs.io/en/latest/)

By  default, both `numpy` and `matplotlib` are installed as dependencies of `obspy`. 

### Installing using pip

You can install `telewavesim` using the [pip package manager](https://pypi.org/project/pip/):

```bash
pip install telewavesim
```
All the dependencies will be automatically installed by `pip`.

### Installing with conda

You can install `telewavesim` using the [conda package manager](https://conda.io).
Its required dependencies can be easily installed with:

```bash
conda install obspy pyfftw -c conda-forge
```

Then `telewavesim` can be installed with `pip`:

```bash
pip install telewavesim
```

#### Conda environment

We recommend creating a custom 
[conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
where `telewavesim` can be installed along with its dependencies. 

- Create a environment called `tws` and install all dependencies:

```bash
conda create -n tws python=3.7 obspy pyfftw -c conda-forge
```

- or create it from the `tws_env.yml` file by first checking out the repository:

```bash
git checkout https://github.com/paudetseis/Telewavesim.git
cd Telewavesim
conda env create -f tws_env.yml
```

Activate the newly created environment:

```bash
conda activate tws
```

Install `telewavesim` with `pip`:

```bash
pip install telewavesim
```

### Installing from source

Download or clone the repository:
```bash
git clone https://github.com/paudetseis/Telewavesim.git
cd Telewavesim
```

Next we recommend following the steps for creating a `conda` environment (see [above](#conda-environment)). Then install using `pip`:

```bash
pip install .
``` 

---
**NOTE**

If you are actively working on the code, or making frequent edits, it is advisable to perform 
installation from source with the `-e` flag: 

```bash
pip install -e .
```

This enables an editable installation, where symbolic links are used rather than straight 
copies. This means that any changes made in the local folders will be reflected in the 
package available on the system.

---

## Usage 

### Jupyter Notebooks

Included in this package is a set of Jupyter Notebooks, which give examples on how to call the various routines and obtain plane wave seismograms and receiver functions. The Notebooks describe how to reproduce published examples of synthetic data from [Audet (2016)](#references) and [Porter et al. (2011)](#references).

- [sim_obs_Audet2016.ipynb](./telewavesim/examples/Notebooks/sim_obs_Audet2016.ipynb): Example plane wave seismograms and P receiver functions for OBS data from [Audet (2016)](#Audet).
- [sim_Prfs_Porter2011.ipynb](./telewavesim/examples/Notebooks/sim_Prfs_Porter2011.ipynb): Example P receiver functions from [Porter et al. (2011)](#Porter)
- [sim_SKS.ipynb](./telewavesim/examples/Notebooks/sim_SKS.ipynb): Example plane wave seismograms for SKS splitting studies.

After [installing `telewaveim`](#installation), these notebooks can be locally installed (i.e., in a local folder `Notebooks`) from the package by running:

```python
from telewavesim import doc
doc.install_doc(path='Notebooks')
```

To run the notebooks you will have to further install `jupyter`:

```bash
conda install jupyter
```

Then ```cd Notebooks``` and type:

```bash
jupyter notebook
```

You can then save the notebooks as `python` scripts, check out the model files and you should be good to go!

### Setting up new models

To set up the models, install the `Jupyter` notebooks and check out the examples in the `models` folder, or visit the [wiki](https://github.com/paudetseis/Telewavesim/wiki/Models) page for `telewavesim`.

### Testing

A series of tests are located in the ``tests`` subdirectory. In order to perform these tests, clone the repository and run `pytest` (`conda install pytest` if needed):

```bash
git checkout https://github.com/paudetseis/Telewavesim.git
cd Telewavesim
pytest -v
```

### API Documentation

The API for all functions in `telewavesim` can be accessed from https://paudetseis.github.io/Telewavesim/.

## References

- Audet, P. (2016). Receiver functions using OBS data: promises and limitations from numerical modelling and examples from the Cascadia Initiative. Geophysical Journal International, 205, 1740-1755. https://doi.org/10.1093/gji/ggw111

- Kennett, B.L.N. (1983). Seismic wave propagation in stratified media. Cambridge University Press, 342pp.

- Porter, R., Zandt, G., & McQuarrie, N. (2011). Pervasive lower-crustal seismic anisotropy in Southern California: Evidence for underplated schists and active tectonics. Lithosphere, 3(3), 201-220. https://doi.org/10.1130/L126.1
