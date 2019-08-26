# Copyright 2019 Pascal Audet
#
# This file is part of Telewavesim.
#
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

> Telewavesim is a Python package for modeling teleseismic (i.e., planar)
body wave propagation through stacks of anisotropic layers for receiver
based studies.

## Licence

# Copyright 2019 Pascal Audet
#
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

## Installation

### Dependencies

You will need **Python 2.7 or 3.5+**.
Also, the following packages are required:

- [`numpy=1.16`](http://numpy.org)
- [`scipy`](https://scipy.org)
- [`pyfftw`](https://github.com/pyFFTW/pyFFTW)
- [`fftw3`](http://www.fftw.org)
- [`obspy`](https://github.com/obspy/obspy/wiki)
__Optional dependencies__ for mapping module and running the Notebooks:

- [`matplotlib`](https://matplotlib.org/)
# - [`pyproj`](https://github.com/jswhit/pyproj)
# - [`cartopy`](https://scitools.org.uk/cartopy/docs/latest/)

### Installing using pip

You can install `telewavesim` using the
[`pip package manager`](https://pypi.org/project/pip/) with either
version of Python:

>>> python2 -m pip install pycurious
>>> python3 -m pip install pycurious

All the dependencies will be automatically installed by `pip`.

### Installing using conda

A more straightforward installation for `pycurious` and all of its
dependencies may be deployed with [Docker](https://www.docker.com).
To install the docker image and start the Jupyter notebook examples:

>>> docker pull brmather/pycurious:latest
>>> docker run --name pycurious -p 8888:8888 brmather/pycurious:latest

## Documenation / Notebooks

Jupyter notebooks that demonstrate the functionality of PyCurious and
some common workflows may be installed to a local directory with

```python
import pycurious
pycurious.install_documentation(path="Notebooks")
```

## Contributing to pycurious

We welcome contributions to `pycurious`, large or small [\\*](#footnote1).
That can be in the form of new code, improvements to the documentation, 
helping with a missing test, or it may even just be pointing out a bug or
potential improvement to the team. 

For bugs and suggestions, the most effective way to reach the team is by
raising an issue on the github issue tracker. Github allows you to classify
your issues so that we know if it is a bug report, feature request or
feedback to the authors.

If you wish to contribute some changes to the code then you should submit a
*pull request*. On GitHub we can review the code that you are contributing
and discuss it before the changes are merged into the development version of
the code. Before embarking on changes to the code, please first take a look
at the [`dev` branch](https://github.com/brmather/pycurious/tree/dev) which
is where unreleased changes are staged: we may have been working on something
similar already. 

### How to create a Pull Request (PR)

Create your own fork of the project on github: log into your github account,
navigate to the [`pycurious` repository](https://github.com/brmather/pycurious/tree/dev)
and click on the *fork* button at the top right corner of this repository.
You can find detailed instructions and a discussion of how to fork a
repository, clone it locally and work on the changes
[in the GitHub guides](https://guides.github.com/activities/forking/). 

### When to make your pull request

It is much easier to merge in small changes to the code than extensive ones
that touch multiple files. If a small incremental change is possible, please
issue a request for that change rather than saving everything up. 

It is a good idea to commit frequently and with commit messages that refer
to the changes that have been make and why they have been made. The commit
message is usually browsed without seeing the actual changes themselves and
which sections of the code have been altered so a helpful message provides
relevant details. 

### Tests

We use [`pytest`](https://pypi.org/project/pytest/) for the unit testing
framework in `pycurious`. In the source directory this means running:

>>> python setup.py test

The existing tests should be passing before you start coding (help us out with
an issue if that is not the case !) and when you have finished. Any new
functionality should also have tests that we can use to verify the code.
It is important that you make it clear if the original tests have had to
change to accomodate new code / functionality.

<a name="footnote1">*</a> _Small changes are our favorites as they are much
easier to quickly merge into the existing code. Proofreading typos, bug fixes,
... the little stuff is especially welcome._

"""

# -*- coding: utf-8 -*-
# from .documentation import install_documentation
# from .grid import CurieGrid, bouligand2009, tanaka1999, maus1995, ComputeTanaka
# from .optimise import CurieOptimise
# from . import mapping
# from . import download
