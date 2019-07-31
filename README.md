pyHiChi
=======

![Language](https://img.shields.io/badge/language-python-orange.svg)
![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20windows-lightgrey.svg)
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)

The project High-Intensity Collisions and Interactions (Hi-Chi) is an open-source collection of Python-controlled tools for performing simulations and data analysis in the research area of strong-field particle and plasma physics. The tools are being developed in C++ and provide high performance using either local or supercomputer resources. The project is intended to offer an environment for testing, benchmarking and aggregative use of individual components, ranging from basic routines to supercomputer codes.

Getting Started
---------------

This repository uses git submodules. To download the code and dependencies please use
```bash
git clone --recursive https://github.com/hi-chi/pyHiChi
```

Alternatively, omit `--recursive` and manually put `pybind11` into the `3rdparty` directory.

The installation process is described below. The installation via `pip-manager` will be implemented in the future.

Installing
----------

### Dependencies

- Python 3: `python` / `python3` and `python-dev` / `python3-dev` available
- [`numba`](https://numba.pydata.org/)
- [`CMake`](https://cmake.org/) 3.1 or higher
- [`pybind11`](https://github.com/pybind/pybind11), comes as a submodule in this repository
- [`fftw3`](http://www.fftw.org/), will be installed automatically in case it is not available
- Additionally to run our examples: [`numpy`](https://numpy.org/) and [`matplotlib`](https://matplotlib.org/)

### On Linux

First, install the [dependencies](###Dependencies). To build the project one needs `gcc` or `icc` supporting C++11.
Run `./build_linux.sh` with the following options:

- `-openmp` to enable OpenMP support (recommended)
- `-fftw` to enable FFTW support
- `-python <path>` to use a non-standard path to Python

After the installation, the binaries will appear in `../bin`. One needs to copy these files to the folder with the Python script to be executed. For example, one can use small tests from the folder `example-tests`. 

### On Windows

For the installation one needs a compiler with C++11 support, for example, Visual Studio 2015 or later.
Run `build_windows.bat` with the following options:

- Windows-only: `/g <generator>` CMake generator name
- `/openmp` to enable OpenMP support (recommended)
- `/fftw` to enable FFTW support
- `/python <path>` to use a non-standard path to Python

After the installation, the binaries will appear in `../bin`. One needs to copy these files to the folder with the Python script to be executed. For example, one can use small tests from the folder `example-tests`. 

Documentation
-------------

To be added soon. Please contact the developers for details meanwhile.

License
-------

pyHiChi is licensed under the **MIT license**. Please refer to [LICENSE](https://github.com/hi-chi/pyHiChi/blob/master/LICENSE).

Authors
-------

### Core Developers and Supervisors

- Dr. Arkady Gonoskov
- [Dr. Iosif Meyerov](https://sites.google.com/site/iosifmeyeroveng/)
- [Valentin Volokitin](https://github.com/ValentinV95)

### Contributions and Thanks

- TBD
