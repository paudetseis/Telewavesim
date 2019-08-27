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

Note that both `numpy` and `matplotlib` are installed as dependencies of `obspy`. See below for full installation details. You also need to download and install the [`fftw3`](http://www.fftw.org) library independently (conda install of `fftw3` is not working with current `gfortran` compiler). 

### Conda environment

We advise creating a custom [conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
where `telewavesim` can be installed along with its dependencies.

Clone the repository:
```bash
git clone https://github.com/paudetseis/Telewavesim.git
cd Telewavesim
```

Create a new environment with all dependencies:
```bash
conda create -n tws python=3.7 obspy pyfftw -c conda-forge
```
or create it from the `tws_env.yml` file:
```bash
conda env create -f tws_env.yml
```
Activate the newly created environment:
```bash
conda activate tws
```

#### Pointing to `fftw3` library

Finally, edit the setup.py file to modify the content of ```extra_link_args``` in the Extension class to point to your compiled `fftw3` library. In the example below the library is install in ```/usr/local/lib```:

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
to perform the pip installation with the ```-e``` flag. This enables an editable installation, where
symbolic links are used rather than straight copies. This means that any changes made in the
local folders will be reflected in the packages available on the system.

## Examples

```bash
cd examples
python sim_Prfs_Porter2011.py
python sim_Prfs_Porter2011_for.py
```