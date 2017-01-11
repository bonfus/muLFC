Local Field Components
======================

LFC is a minimal library to facilitate the evaluation of the local field at the muon site. 
The original code comes from the muesr repository.


Installing LFC
--------------
This section provides an overview and guidance for installing LFC on
various target platforms.

### Prerequisites

The LFC C library have no dependencies. However, cmake is
needed to compile the code.

The python extension has the following dependencies:

* Python 2.7, 3.1+      (http://www.python.org)
* Numpy 1.6.0+          (http://www.numpy.org)

(Windows is likely to be supported with minor changes to the code.)


### Installation

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
