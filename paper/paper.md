---
title: 'Telewavesim: Python software for teleseismic body wave modeling'
tags:
  - Python
  - Fortran
  - geophysics
  - seismology
  - modeling
authors:
  - name: Pascal Audet
    orcid: 0000-0003-2364-9454
    affiliation: 1
  - Colin J. Thomson
    affiliation: "2, 3"
  - Michael G. Bostock
    orcid: 0000-0003-1172-7240
    affiliation: 4
affiliations:
  - name: Department of Earth and Environmental Sciences, University of Ottawa, Canada
    index: 1
  - name: Sherlock Seismic Solutions, Cambridge, U.K.
    index: 2
  - name: Formerly Queen's University, Canada, and Schlumberger Cambridge Research, U.K.
    index: 3
  - name: Department of Earth, Ocean and Atmospheric Sciences, The University of British Columbia, Canada
    index: 4


date: 2 October 2019
bibliography: paper.bib
---

# Summary

The structure of the Earth's crust and upper mantle provides useful information on the 
internal composition and dynamics of our planet. Some of the most widely used techniques
to infer these properties are based on examining the effect of teleseismic body wave 
(i.e., P and S waves that originate from distant earthquakes and arrive as plane waves)
propagation through stratified media. Modeling the 
seismic response from stacks of subsurface layers is therefore an essential tool in 
characterizing their effect on observed seismograms. Although a few established techniques
are available for this purpose (e.g., @Frederiksen:2000), these are either based on 
approximations or else do not handle general anisotropy or the effect of an overlying
fluid medium.

``telewavesim`` is a Python package for synthesizing teleseismic 
body-wave propagation through stacks of generally anisotropic and strictly horizontal 
layers using the matrix propagator approach of @Kennett:1983 and implemented by
@Thomson:1997. Python 
allows wrapping low-level Fortran routines for speed while maintaining flexibility
and ease of use for code interaction and plotting. The integration of ``obspy`` 
[@obspy] ``Stream`` objects allows flexibility in manipulating and visualzing 
the resulting seismograms. The software contains a wide range of elastic stiffness 
definitions extracted from the literature (@Brownlee:2017) to accurately represent 
seismic anisotropy due to various mineralogies and rock types. The software also 
accurately models acoustic reverberations from an overlying column of water, 
effectively simulating ocean-bottom seismograph (OBS) station recordings. 

The software will be useful in a variety of teleseismic receiver-based studies, 
such as P or S receiver functions, shear-wave splitting from 
core-refracted teleseismic shear waves (i.e., SKS, SKKS), etc. It may also be the 
starting point for stochastic inverse methods (e.g., Monte Carlo sampling) and more
general (i.e., point) sources through slowness integration (@Frazer:1984}. Common 
computational workflows that reproduce published examples [@Audet:2016; @Porter:2011] 
are covered in the Jupyter notebooks bundled with this package. These notebooks 
can be further used in a teaching environment to study the effects of seismic wave 
scattering effects in stratified media. The API documentation is up-to-date on 
[GitHub pages](https://paudetseis.github.io/Telewavesim/). `telewavesim` and all Python
dependencies can be installed through the pypi.org `pip` package manager.

# Acknowledgements

This work was supported by the Natural Sciences and Engineering Research Council of Canada.

# References