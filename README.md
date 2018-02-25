Local Field Components
======================

LFC is a small library to facilitate the evaluation of the local field at the muon site. 


Install
-------

The LFC C library has no dependencies. CMake is used as a building tool.
Any C90 compiler can be used to build the library.
The details for the Python package follows.

### Prerequisites

The python extension has the following dependencies:

* Python 2.7, 3.1+      (http://www.python.org)
* Numpy 1.8.0+          (http://www.numpy.org)


### Installation

The easiest way to install the Python extension is using pre-baked
wheels with pip:

    pip install mulfc

**Important note:** when installing on *Windows* with `pip`, the minimal
dependencies for numpy are:

| Python version | Numpy version |
|----------------|---------------|
| 2.7            | 1.8+          |
| 3.5            | 1.10+         |
| 3.6            | 1.12+         |

Python versions are different from the one listed above are not
distributed with precompiled wheels.

### Compilation

You may want to build the extension yourself, especially for optimizing performances.

To compile and install the C library just do

```bash
mkdir build
cd build
cmake ..
make
make install
```

To install the python extension you can do

```bash
cd python
python setup.py build
python setup.py test
python setup.py install
```

The library and the python extension are independent from each other.

Usage
-----

See documentation of the various functions.
