#!/usr/bin/python
"""
 Usage: [input file]

 Below we provide examples of custom usage
>>> 
>>> from Coulomb_t import Coulomb
>>> run = Coulomb()
>>> task = Task(task='camm', trans=False, mol=<Molecule>, basis_set='6-31G*', **kwargs)
>>> run.eval(task)
>>> run.save(property='dma', out='mol.par')
>>> print run
>>> run.clear()
>>> 
"""
# ---------------------------
from sys import argv, exit
print __doc__
if len(argv)==1: exit(0)
# ---------------------------------
from Coulomb_t import Coulomb, Task
# ---------------------------------
def parse_input(file): 
    """Parser for Coulomb input files. This function is supposed to return a list of Task instances"""
    test_task = Task()
    return [test_task]*3

# instantiate the Coulomb routine
run   = Coulomb()

# parse the input file
tasks = parse_input(argv[1])

# run Coulomb
for task in tasks:
    run.eval(task)

