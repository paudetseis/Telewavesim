cmake_minimum_required(VERSION 3.17.2...3.29)
project(${SKBUILD_PROJECT_NAME} LANGUAGES Fortran)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

# F2PY headers
execute_process(
  COMMAND "${PYTHON_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")

# Find the FFTW3 installation on your system
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW)
link_libraries(PkgConfig::FFTW)

target_link_libraries(fortranobject PUBLIC  Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

add_custom_command(
  OUTPUT mccoremodule.c mccore-f2pywrappers.f
  DEPENDS src/relmt/mccore.f90
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          "${CMAKE_CURRENT_SOURCE_DIR}/src/relmt/mccore.f90" -m mccore --lower)

python_add_library(
  mccore MODULE "${CMAKE_CURRENT_BINARY_DIR}/mccoremodule.c"
  "${CMAKE_CURRENT_BINARY_DIR}/mccore-f2pywrappers.f"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/relmt/mccore.f90" WITH_SOABI)
target_link_libraries(mccore PRIVATE fortranobject)

install(TARGETS mccore DESTINATION .)