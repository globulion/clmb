clmb
====

Coulomb.py Project

The project was lauched in 2010 in Wrocław University of Science and Technology
in Wrocław, Poland.
The package is designed for quantum chemistry calculations
of electrostatic interaction energy between two molecular species.
Currently, the implemented methods are:

* Electrostatic energy from distributed charges,
* Electrostatic energy from distributed multipoles,
* Electrostatic energy from density cube files,
* Electrostatic energy from exact first-order expression for two interacting charge densities. 
   
The auxiliary methods implemented are:

* `ESP` - charges from fitting to the electrostatic potential,
* `CAMM` - Cumulative Atomic Multipole Moments.
* `CABMM` - Cumulative Atomic and Bond Moments

The tutorial is under preparation. 

### Installation ###

Installation prerequisites:
- Python 2.6 or newer
- NumPy module for Python
- pp module for Python
- PyQuante modified version
- LibInt 1.1.4

To install the Coulomb package type the following commands:
```
sudo python setup.py install
```
