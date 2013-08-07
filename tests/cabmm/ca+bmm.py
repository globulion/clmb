#!/usr/bin/python

import coulomb.work
from coulomb import *
from coulomb.run import *
import sys

h2o = Molecule('H2O',
                   [(8,  ( 0.00000000,     0.00000000,     0.04851804)),
                    (1,  ( 0.75300223,     0.00000000,    -0.51923377)),
                    (1,  (-0.75300223,     0.00000000,    -0.51923377))],
                   units='Angstrom')

#bonds = [ (1,0), (2,0) ]
#bonds = [(4,5)]
bonds = []
# calculate camms
basis=sys.argv[1]
c = multip.MULTIP(h2o,basis,'HF',bonds=bonds)
c.camms()
c.__printCAMMs__()
print c.Mon
r = zeros(3)
for i in c.Dip: r+=i
print r
print sum(c.Mon)
c.mmms()
c.__printMMMs__()
D = c.get()
print c
#for j in D: print j.OverallMoments()
