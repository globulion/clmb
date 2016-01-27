#!/usr/bin/python
"""
 Tests (just sketch - not run yet)
"""
import os, filecmp
from Coulomb_t import Coulomb, Task

# file with xyz structure
xyz_file = '%s/data/water.xyz' % os.environ['PWD']

# load Coulomb
C = Coulomb()

# evaluate
C('camm', molecules=xyz_file, basis='6_31G*').eval()

# save the output
C.dma.write('water.camm')


# --- alternative way
# load task
T = Task('camm', molecules=xyz_file, basis='6_31G*')
C.eval(T)
C.dma.write('water.camm2')

# check if the output is the same
assert filecmp.cmp('water.camm', 'water.camm2')
os.system('rm water.camm2')
