project('telewavesim', 'c',
  version : '0.3.0',
  license: 'MIT',
  meson_version: '>=0.64.0',
  default_options : [
    'warning_level=1',
    'c_args=-Wno-unused-function -Wno-conversion -Wno-misleading-indentation -Wno-incompatible-pointer-types',
    'fortran_args=-Wno-conversion',
    'fortran_std=legacy'],
)

add_languages('fortran', native: false)

py_mod = import('python')
py = py_mod.find_installation()
py_dep = py.dependency()

lapack_dep = dependency('lapack', required: false)

incdir_numpy = run_command(py,
  ['-c', 'import os; os.chdir(".."); import numpy; print(numpy.get_include())'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; os.chdir(".."); import numpy.f2py; print(numpy.f2py.get_include())'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py)

# Unlike distutils, meson doesn't yet include some of the f2py stuff
fortranobject_c = incdir_f2py / 'fortranobject.c'

fortranobject_lib = static_library('_fortranobject',
  fortranobject_c,
  dependencies: py_dep,
  include_directories: [incdir_numpy, incdir_f2py])

fortranobject_dep = declare_dependency(
  link_with: fortranobject_lib,
  include_directories: [incdir_numpy, incdir_f2py])


rmat_source = custom_target('rmat_fmodule.c',
  input : ['src/rmat_sub.f90', 'src/rmat.f90'], 
  output : ['rmat_fmodule.c', 'rmat_f-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py', '@INPUT@', '-m', 'rmat_f', '--lower']
)

py.extension_module('rmat_f',
  [
    'src/rmat_sub.f90',
    'src/rmat.f90', 
    rmat_source,
  ],
  include_directories: inc_np, 
  dependencies : [py_dep, fortranobject_dep, lapack_dep],
  subdir: 'telewavesim',
  install : true
)

py.install_sources(
  'telewavesim/utils.py',
  'telewavesim/elast.py',
  'telewavesim/wiggle.py',
  'telewavesim/doc.py',
  'telewavesim/__init__.py',
  pure: false,
  subdir: 'telewavesim'
  )

install_data(
    'telewavesim/examples/models/model_Audet2016.txt',
    'telewavesim/examples/models/model_Porter2011.txt',
    'telewavesim/examples/models/model_SKS.txt',
    install_dir: py.get_install_dir(pure: false) / 'telewavesim/examples/models',
)
install_data(
    'telewavesim/examples/Notebooks/sim_obs_Audet2016.ipynb',
    'telewavesim/examples/Notebooks/sim_Prfs_Porter2011.ipynb',
    'telewavesim/examples/Notebooks/sim_SKS.ipynb',
    install_dir: py.get_install_dir(pure: false) / 'telewavesim/examples/Notebooks',

)