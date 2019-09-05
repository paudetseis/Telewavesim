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

The current version was developed using **Python3.7** \
Also, the following packages are required:

- ``gfortran`` (https://gcc.gnu.org/wiki/GFortran) (or any Fortran compiler)
- ``obspy`` (https://github.com/obspy/obspy/wiki)
- ``pyfftw`` (https://pyfftw.readthedocs.io/en/latest/)

By  default, both ``numpy`` and ``matplotlib`` are installed as dependencies of ``obspy``. \
See below for full installation details. 

Download the software
+++++++++++++++++++++

- Clone the repository:

.. sourcecode:: bash

   git clone https://github.com/paudetseis/Telewavesim.git
   cd Telewavesim

Conda environment
+++++++++++++++++

We recommend creating a custom ``conda`` environment
where ``telewavesim`` can be installed along with its dependencies.

.. sourcecode:: bash

   conda create -n tws python=3.7 obspy pyfftw -c conda-forge

or create it from the ``tws_env.yml`` file:

.. sourcecode:: bash

   conda env create -f tws_env.yml

Activate the newly created environment:

.. sourcecode:: bash

   conda activate tws

Installing using pip
++++++++++++++++++++

Once the previous steps are performed, you can install ``telewavesim`` using pip:

.. sourcecode:: bash

   pip install .

.. note::

   Please note, if you are actively working on the code, or making frequent edits, it is advisable
   to perform the pip installation with the ``-e`` flag. This enables an editable installation, where
   symbolic links are used rather than straight copies. This means that any changes made in the
   local folders will be reflected in the packages available on the system.


"""
# -*- coding: utf-8 -*-
# from .doc import install_doc
# from . import elast
# from . import conf
# from . import utils
# from . import wiggle
# from .rmat_f import *
