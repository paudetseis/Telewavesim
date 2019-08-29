# Copyright 2019 Pascal Audet

# This file is part of Telewavesim.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""

Telewavesim is a Python package for modeling teleseismic (i.e., planar)
body wave propagation through stacks of anisotropic layers for receiver
based studies.

Licence
-------

Copyright 2019 Pascal Audet

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Installation
------------

Dependencies
++++++++++++

The current version was developed using **Python3.7** \
Also, the following packages are required:

- ``obspy`` (https://github.com/obspy/obspy/wiki)
- ``pyfftw`` (https://pyfftw.readthedocs.io/en/latest/)
- ``fftw`` (http://www.fftw.org)

By  default, both ``numpy`` and ``matplotlib`` are installed as dependencies of ``obspy``. \
See below for full installation details. 

Conda installation
++++++++++++++++++

We advise creating a custom ``conda`` environment
where ``telewavesim`` can be installed along with its dependencies.

Clone the repository:

.. sourcecode:: bash

   git clone https://github.com/paudetseis/Telewavesim.git
   cd Telewavesim

Create a new environment and install all dependencies:

.. sourcecode:: bash

   conda create -n tws python=3.7 obspy pyfftw -c conda-forge

or create it from the ``tws_env.yml`` file:

.. sourcecode:: bash

   conda env create -f tws_env.yml

Activate the newly created environment:

.. sourcecode:: bash

   conda activate tws

Fortran compilation and ``fftw3`` library
*****************************************

You can further use ``conda`` to install the required Fortran compiler and the \
``fftw`` library. This is the default installation, as the ``os.environ`` points \
to the ``tws`` environment library for dynamic linking. In this case, install \
``gfortran`` and ``fftw`` using ``conda``. On a MacOSX, the ``gfortran`` package \
is ``gfortran_osx-64``; for Linux, the ``gfortran`` package is ``gfortran_linux-64`` \
(check out https://anaconda.org/search?q=gfortran for the available packages):

.. sourcecode:: bash

   conda install gfortran_osx-64 fftw

You can check that the active Fortran compiler resides in the ``tws`` environment:

.. sourcecode:: bash

   which gfortran


Separate Fortran build
**********************

If you wish to use a different Fortran compiler available system-wide (e.g., \
Intel Fortran installed in ``/usr/local/bin``), you will need to independently \
download and install the ``fftw`` library (http://www.fftw.org), then edit the 
``setup.py`` file to modify the content of ``extra_link_args`` in the ``Extension`` \
class to point to the location of your compiled ``fftw3`` library. In the example 
below the library is installed in ``/usr/local/lib``:

.. sourcecode:: python

   ext = [Extension(name='telewavesim.rmat_f',
                    sources=['src/rmat.f90', 'src/rmat_sub.f90'],
                    extra_f90_compile_args=["-O3"],
                    extra_link_args=["-L/usr/local/lib", "-lfftw3"])]


1) Installing using pip
+++++++++++++++++++++++

Once the previous steps are performed, you can install ``telewavesim`` using pip:

.. sourcecode:: bash

   pip install .

2) Building and Installing
++++++++++++++++++++++++++

Alternatively, you can build and install the project (from the root of the source tree, e.g., inside the cloned `telewavesim` directory):

.. sourcecode:: bash

   python setup.py build 
   python setup.py install


Please note, if you are actively working on the code, or making frequent edits, it is advisable
to perform the pip installation with the ``-e`` flag. This enables an editable installation, where
symbolic links are used rather than straight copies. This means that any changes made in the
local folders will be reflected in the packages available on the system.


"""