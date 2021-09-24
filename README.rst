pyHiChi
=======


.. |Language| image:: https://img.shields.io/badge/language-python-orange.svg
.. |Platform| image:: https://img.shields.io/badge/platform-linux%20%7C%20windows-lightgrey.svg
.. |License| image:: https://img.shields.io/badge/license-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT

|Language| |Platform| |License|


The project High-Intensity Collisions and Interactions (Hi-Chi) is an open-source collection of Python-controlled tools for performing simulations and data analysis in the research area of strong-field particle and plasma physics. The tools are being developed in C++ and provide high performance using either local or supercomputer resources. The project is intended to offer an environment for testing, benchmarking and aggregative use of individual components, ranging from basic routines to supercomputer codes.


Getting Started
---------------

This repository uses git submodules. To download the code and dependencies please use

.. code:: bash

  git clone --recursive https://github.com/hi-chi/pyHiChi


Alternatively, omit ``--recursive`` and manually put ``pybind11`` into the ``3rdparty`` directory.

The installation process is described below. The installation via ``pip-manager`` will be implemented in the future.


Installing
----------

.. _dependencies:

Dependencies
^^^^^^^^^^^^
.. reST doesn't allow nested markup. The following is a workaround that allows the link texts to be effectively nested with a modifier, using "roles".
.. _numba: https://numba.pydata.org/
.. _cmake: https://cmake.org/
.. _pybind11: https://github.com/pybind/pybind11
.. _fftw3: http://www.fftw.org/
.. _numpy: https://numpy.org/
.. _matplotlib: https://matplotlib.org/

.. |numba| replace:: ``numba``
.. |cmake| replace:: ``CMake``
.. |pybind11| replace:: ``pybind11``
.. |fftw3| replace:: ``fftw3``
.. |numpy| replace:: ``numpy``
.. |matplotlib| replace:: ``matplotlib``


- Python 3: ``python`` / ``python3`` and ``python-dev`` / ``python3-dev`` available
- |numba|_
- |cmake|_ 3.1 or higher.
- |pybind11|_, comes as a submodule in this repository.
- |fftw3|_, will be installed automatically in case it is not available.
- Additionally to run our examples: |numpy|_ and |matplotlib|_.

The easiest way to install the required Python libraries is to run

.. code:: bash

  pip install -r requirements.txt

The ``requirements.txt`` file contains a complete list of the necessary Python libraries, including version control where needed.


On Linux/Mac
^^^^^^^^^^^^
First, install the :ref:`dependencies`. To build the project one needs ``gcc`` or ``icc`` supporting C++11.
Run ``./build_linux.sh`` with the following options:

- ``-openmp`` to enable OpenMP support (recommended)
- ``-fftw`` to enable FFTW support
- ``-python <path>`` to use a non-standard path to Python

After the installation, the binaries will appear in ``../bin``. One needs to copy these files to the folder with the Python script to be executed. For example, one can use small tests from the folder ``example-tests``.


On Windows
^^^^^^^^^^
For the installation one needs a compiler with C++11 support, for example, Visual Studio 2015 or later.
Run ``build_windows.bat`` with the following options:

- Windows-only: ``/g <generator>`` CMake generator name
- ``/openmp`` to enable OpenMP support (recommended)
- ``/fftw`` to enable FFTW support
- ``/python <path>`` to use a non-standard path to Python

After the installation, the binaries will appear in ``../bin``. One needs to copy these files to the folder with the Python script to be executed. For example, one can use small tests from the folder ``example-tests``.


Documentation
-------------

To be added soon. Please contact the developers for details meanwhile.


License
-------

pyHiChi is licensed under the **MIT license**. Please refer to `LICENSE <https://github.com/hi-chi/pyHiChi/blob/master/LICENSE>`_.


Authors
-------


Core Developers and Supervisors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

- Dr. Arkady Gonoskov
- `Dr. Iosif Meyerov <https://sites.google.com/site/iosifmeyeroveng/>`_
- `Valentin Volokitin <https://github.com/ValentinV95>`_


Contributions and Thanks
^^^^^^^^^^^^^^^^^^^^^^^^

- TBD
