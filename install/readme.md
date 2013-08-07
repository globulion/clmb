Installation steps
==================

Coulomb.py requires many prerequisites:
* Python 2.7.3 or newer but not Python 3xxx
* SciPy 0.9.0 or newer
* setuptools module for Python
* python-dev
* libbbg (download it [here](https://github.com/globulion/libbbg))
* PyQuante Modified Version (see this directory)

Also, recommended are:
* libint 1.1.4 (see this directory)

The procedure for installation is as follows:
* install Python, Scipy, scitools, setuptools, python-dev and [libbbg](https://github.com/globulion/libbbg "Click here to get the info how to install libbbg")
* install libint (if want to increase the speed of ERI computation)
* install PyQuante Modified Version
* install coulomb

## I. Installation of libint 1.1.4


Three steps are generally necessary:
```
./configure --enable-shared
make
sudo make install
```
If something is not working (it should be ok) change the settings in `configure`
and make sure you have installed `g++` compiler. The installation takes a longer while.

## II. Installation of PyQuante Modified Version


First, you have to set appropriate paths in file `setup.py` (for libint directory
if it was installed).  
Then run the following command:
```
sudo python setup.py install --enable-libint
```
Then probably a reboot of the system is necessary. Try after reboot to run Python
and import PyQuante module:
``` import PyQuante```.
If everything is ok no message is printed. You can however get two messages:
* module `cints` is not installed
* openbabel is not installed (it is not harmful)

The first problem is due to the failure in compiling and linking extension module `cints`. In this case check
the compiler and so on.

## III. Installation of Coulomb.py

Just type the following command:
```
sudo python setup.py install
```

Everything should work properly. Now run the test computation (you should be in coulomb main directory):
```
cd tests/ethylene
clmb hf.inp
clmb e2-eds.inp
clmb zwykle_e2.inp
```
If everything is OK you can be happy now.
