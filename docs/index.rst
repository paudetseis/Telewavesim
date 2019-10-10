
.. figure:: ../telewavesim/examples/picture/tws_logo.png
   :align: center

Telewavesim documentation
=========================

This package contains ``python`` and ``fortran`` modules to synthesize teleseismic body-wave propagation through stacks of generally anisotropic and strictly horizontal layers using the matrix propagator approach of `Kennett (1983) <https://www.oapen.org/search?identifier=459524>`_, as implemented in `Thomson (1997) <https://doi.org/10.1016/S0031-9201(97)00033-2>`_. The software also properly models reverberations from an overlying column of water using the R/T matrix expressions of `Bostock and Trehu (2012) <https://doi.org/10.1785/0120110162>`_, effectively simulating ocean-bottom seismic (OBS) station recordings. The software will be useful in a variety of teleseismic receiver-based studies, such as P or S receiver functions, long-period P-wave polarization, shear-wave splitting from core-refracted shear waves (i.e., SKS, SKKS), etc. It may also be the starting point for stochastic inverse methods (e.g., Monte Carlo sampling). The main part of the code is written in ``fortran`` with ``python`` wrappers. Common computational workflows are covered in the ``Jupyter`` notebooks bundled with this package.

Quick links
"""""""""""

* `Git repositories <https://github.com/paudetseis/Telewavesim>`_

.. toctree::
   :maxdepth: 1
   :caption: Getting Started

   api

.. toctree::
   :maxdepth: 1
   :caption: Module

   elast
   utils
   wiggles

.. toctree::
   :maxdepth: 1
   :caption: Jupyter Notebooks

   Example 1: A simple example for OBS <https://github.com/paudetseis/Telewavesim/blob/master/telewavesim/examples/Notebooks/sim_obs_Audet2016.ipynb>
   Example 2: Lower crustal anisotropy <https://github.com/paudetseis/Telewavesim/blob/master/telewavesim/examples/Notebooks/sim_Prfs_Porter2011.ipynb>
   Example 3: Simulating SKS splitting <https://github.com/paudetseis/Telewavesim/blob/master/telewavesim/examples/Notebooks/sim_SKS.ipynb>

