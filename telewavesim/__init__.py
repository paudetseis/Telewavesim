# Copyright 2019 Pascal Audet

# This file is part of Telewavesim.

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

"""

Telewavesim is a software for modeling teleseismic (i.e., planar)
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

The current version works well with **Python>3.5**.
Also, the following packages are required:

- `gfortran <https://gcc.gnu.org/wiki/GFortran>`_ (or any Fortran compiler)
- `obspy <https://github.com/obspy/obspy/wiki>`_

By  default, both ``numpy`` and ``matplotlib`` are installed as dependencies
of ``obspy``.

Conda environment
+++++++++++++++++

We recommend creating a custom
`conda environment <https://conda.io/docs/user-guide/tasks/manage-environments.html>`_
where ``telewavesim`` can be installed along with its dependencies:

.. sourcecode:: bash

   conda create -n tws python=3.7 obspy -c conda-forge

Activate the newly created environment:

.. sourcecode:: bash

   conda activate tws

Installing latest version from PyPi
+++++++++++++++++++++++++++++++++++

Install the latest version from PyPi with the following command:

.. sourcecode:: bash

    pip install telewavesim


Installing development version from source
++++++++++++++++++++++++++++++++++++++++++

- Clone the repository:

.. sourcecode:: bash

   git clone https://github.com/paudetseis/Telewavesim.git
   cd Telewavesim

- Install using ``pip``:

.. sourcecode:: bash

   pip install .

Possible installation pitfalls with conda
+++++++++++++++++++++++++++++++++++++++++

Using ``conda`` it might be necessary to use the fortran compiler provided with
conda-forge. Add ``gfortran_osx-64`` or ``gfortran_linux-64`` package to the
above ``conda`` environment calls.
On Linux it might further be necessary to install the ``lapack`` conda package.

Usage
-----

Jupyter Notebooks
+++++++++++++++++

Included in this package is a set of Jupyter Notebooks (see Table of Content),
which give examples on how to call the various routines and obtain plane wave
seismograms and receiver functions.
The Notebooks describe how to reproduce published examples of synthetic data
from `Audet (2016) <https://doi.org/10.1093/gji/ggw111>`_ and
`Porter et al. (2011) <https://doi.org/10.1130/L126.1>`_.


After ``telewaveim``, these notebooks can be locally installed
(i.e., in a local folder ``Notebooks``) from the package
by typing in a ``python`` window:

.. sourcecode:: python

   from telewavesim import doc
   doc.install_doc(path='Notebooks')

To run the notebooks you will have to further install ``jupyter``.
From the terminal, type:

.. sourcecode:: bash

   conda install jupyter

Followed by:

.. sourcecode:: bash

   cd Notebooks
   jupyter notebook

You can then save the notebooks as ``python`` scripts,
check out the model files and you should be good to go!

Setting up new models
+++++++++++++++++++++

To set up the models, install the ``Jupyter`` notebooks and
check out the examples in the ``models`` folder, or visit the
`wiki <https://github.com/paudetseis/Telewavesim/wiki/Models>`_ page
for ``telewavesim``.

Testing
+++++++

A series of tests are located in the ``tests`` subdirectory.
In order to perform these tests, run ``pytest``
(``conda install pytest`` if needed):

.. sourcecode:: bash

   pytest -v --pyargs telewavesim



"""

__version__ = '0.2.0-dev'
