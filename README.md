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
* Numpy 1.6.0+          (http://www.numpy.org)


### Installation

The easiest way to install the python extension muLFC is using pre-baked
wheels with pip:

    pip install mulfc
    
or packagecloud:

[packagecloud.io/muLFC/wheels](https://packagecloud.io/muLFC/wheels)

In this latter case you have to choose the right package for your python installation. 
For example, `LFC-0.1-cp36-cp36m-manylinux1_i686.whl` is compatible with
Python 3.6 (cp36) and a 32 bit system, while `LFC-0.1-cp34-cp34m-manylinux1_x86_64.whl`
is for Python 3.4 installed on a 64 bit system.

### Compilation

You may want to build the code yourself, especially for optimizing performances.

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
