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
- ``fftw3`` (http://www.fftw.org)

By  default, both ``numpy`` and ``matplotlib`` are installed as dependencies of ``obspy``. \
See below for full installation details. You also need to download and install the \
``fftw3`` library independently (conda install of ``fftw3`` is not currently working with \
independent GCC's ``gfortran`` build). 

Conda environment
+++++++++++++++++

We advise creating a custom [conda environment](https://conda.io/docs/user-guide/tasks/manage-environments.html)
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

Pointing to ``fftw3`` library
***************************

Finally, edit the `setup.py` file to modify the content of ``extra_link_args`` in the Extension class to point to your compiled ``fftw3`` library. In the example below the library is install in ``/usr/local/lib``:

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