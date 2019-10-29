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
Use the `telewaveim.doc.install_doc` function to copy all
Jupyter Notebooks and example data to a local directory.

"""

import pkg_resources as _pkg_resources
from distutils import dir_util as _dir_util


def install_doc(path="./Telewavesim-Examples"):
    """
    Install the examples for Telewavesim in the given location.

    WARNING: If the path exists, the files will be written into the path
    and will overwrite any existing files with which they collide. The default
    path ("./Telewavesim-Examples") is chosen to make collision less
    likely/problematic.

    The documentation for Telewavesim is in the form of jupyter notebooks.

    """

    Notebooks_Path = _pkg_resources.resource_filename(
        "telewavesim", "examples")

    _dir_util.copy_tree(
        Notebooks_Path,
        path,
        preserve_mode=1,
        preserve_times=1,
        preserve_symlinks=1,
        update=0,
        verbose=1,
        dry_run=0,
    )
