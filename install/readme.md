Installation steps
==================

Coulomb.py requires many prerequisites:
* SciPy 0.9.0 or newer
* libbbg (download it [here](https://github.com/globulion/libbbg))
* PyQuante Modified Version (see this directory)

Also, recommended are:
* libint 1.1.4 (see this directory)

The procedure for installation is as follows:
* install Scipy, scitools and [libbbg](https://github.com/globulion/libbbg,"r")
* install libint
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

