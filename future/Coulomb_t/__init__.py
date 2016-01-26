
__version__ = "0.0.0"  # NOTE: keep in sync with version in setup.py

__all__ = ['Coulomb', 'Task'] 

class Task(object):
   """Contains the specification of the task to perform by Coulomb"""

   def __init__(self, name=None, molecules=list(), transition=None, **kwargs): 
       self._name      = name
       self._molecules = molecules
       # support passing single molecule without the need of putting it inside a list
       if not isinstance(molecules, list): self._molecules = list([self._molecules])  
       self._transition= transition
       # set the method-specific options
       self._set_opts(name, **kwargs)

   def _set_opts(self, name, **kwargs):
       """Set the options according to the method"""
       if   name == 'camm'  : pass
       elif name == 'esp'   : pass
       elif name == 'eleds' : pass
          
   # this sets the molecule or displays the current list of added molecules
   @property
   def molecule(self): 
       return self._molecules
   @molecule.setter
   def molecule(self, mol):
       self._check_molecule(mol)
       self._molecules.append(mol)
       return

   @property
   def name(self): return self._name
   @name.setter
   def name(self, name): self._name = name

   # helper methods

   def _check_molecule(self, mol):
       """Checks if the molecule is fully specified"""
       assert mol.charge       is not None, " Error: Molecule %s has not been assigned its charge!      " % mol.name
       assert mol.multiplicity is not None, " Error: Molecule %s has not been assigned its multiplicity!" % mol.name
       assert mol.basis_set    is not None, " Error: Molecule %s has not been assigned its basis set!   " % mol.name
       return    

   def __repr__(self):
       """Print me!"""
       log = 'A happy task!'
       return str(log)

# ------------------------------
from core.Coulomb import Coulomb
# ------------------------------

#  info: 
#  core.Coulomb is a main module (driving head)
#  core.Coulomb::Coulomb is a class that contains the driving routines

