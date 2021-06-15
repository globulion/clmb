# ----------------------------------------------------------- #
#               RUN - BASIC OPERATION PERFORMER               #
# ----------------------------------------------------------- #

#---------------------------------
from numpy import *
from numpy.linalg import svd
#---------------------------------
import re, os, sys, getopt, random
from time import time
from sys  import argv, exit
#------------------------------------------------------------------------
from PyQuante.Molecule import Molecule
from PyQuante          import SCF
from PyQuante.Ints     import getbasis, getS, getM, getJ, getK, sortints
from PyQuante.cints    import ijkl2intindexBTF as intindexBTF
from PyQuante.CGBF     import coulomb
#------------------------------------------------------------------------
from libbbg.utilities import *
from util import *

class RUN:
  """calculate basic variables and matrices""" 
  
  def __init__(self,molecule,basis,method,matrix=0,multInts=0,hexadecapoles=False):
     # start measuring time of action
     self.clock = TIMER()
     # molecule
     self.molecule = molecule 
     self.N_at = len(molecule.atoms)
     self.RArray = zeros((self.N_at,3),dtype=float64)
     for i in range(self.N_at):
         self.RArray[i,:] = array(self.molecule.atoms[i].pos())
     # basis set
     self.bfs = getbasis(molecule,basis)
     self.K = len(self.bfs)
     self.basis = basis
     self.method = method
     self.matrix = matrix
     self.multInts = multInts

     if ('NoneType' in str(type(matrix)) and 'NoneType' in str(type(multInts))):
        ### 1-enectron integrals and 2-electron integrals
        self.solv= SCF(molecule,method="%s"%method,basis="%s"%basis)
        self.clock.actualize('one- and two-electron integrals')
        ### start SCF calculations
        self.solv.iterate()
        self.clock.actualize('SCF calculations')

        # energy    
        self.E = self.solv.energy
        # density matrix
        self.P = 2*self.solv.dmat
        # overlap matrix
        self.S = self.solv.S
     elif ('array' in str(type(matrix)) and 'NoneType' in str(type(multInts))):
        self.P = matrix#matrix['P']
        #self.solv= SCF(molecule,method="%s"%method,basis="%s"%basis)
        #self.solv.iterate()
        #self.S = self.solv.S
        self.S = getS(self.bfs)
     elif ('array' in str(type(matrix)) and 'array' in str(type(multInts))):
        self.P = matrix#matrix['P']
    
     # calculate multipole integrals 
     if not multInts:
        self.D,self.Q,self.O,self.H = getM(self.bfs,hexadecapoles=hexadecapoles) 
        self.clock.actualize('multipole integrals evaluation')
     else: 
        self.S,self.D,self.Q,self.O,self.H = multInts 
        
        
   
   

### 
### The 'run' class below is during preparation. It was thought to be
### a singleton class for performing multiple tasks (e.g. interaction
### energy) for the same pair of molecules but running basic operands
### as density matrix and multipole integrals only once. This is not
### the case for actual implementation of COULOMB when the program
### calculates these working variables each time for each task speci-
### fied in the input.
### 

# ---------------------------------------------- #
#      RUN - A BASIC SINGLETON FOR EACH TASK     #
# -----------------------------------------------#

class START:
  """Singleton class"""
  
  def __init__(self): pass
  a = 0
  @staticmethod
  def start():
   if not START.a:
     print("Wykonuję się tylko JEDEN raz i ustawiam a na 1!")
     START.a = 1
     
class run:
  """calculate basic variables and matrices""" 
  
  def __init__(self): pass
  start = False
  @staticmethod
  def run(molecule,basis,method,matrix=None,multInts=None):
   if not run.start: 
     print("Calculating basic variables")
     # start measuring time of action
     run.clock = TIMER()
     # molecule
     run.molecule = molecule
     run.N_at = len(molecule.atoms)
     run.RArray = zeros((run.N_at,3),dtype=float64)
     for i in range(run.N_at):
         run.RArray[i,:] = array(run.molecule.atoms[i].pos())
     # basis set
     run.bfs = getbasis(molecule,basis)
     run.K = len(run.bfs)
     run.basis = basis
     run.method = method
     run.matrix = matrix
     run.multInts = multInts

     if (not matrix and not multInts):
        ### 1-enectron integrals and 2-electron integrals
        run.solv= SCF(molecule,method="%s"%method,basis="%s"%basis)
        run.clock.actualize('one- and two-electron integrals')
        ### start SCF calculations
        run.solv.iterate()
        run.clock.actualize('SCF calculations')

        # energy    
        run.E = run.solv.energy
        # density matrix (RHF)
        run.P = 2*run.solv.dmat
        # overlap matrix
        run.S = run.solv.S
     elif (matrix and not multInts):
        run.P = matrix['P']
        run.solv= SCF(molecule,method="%s"%method,basis="%s"%basis)
        run.solv.iterate()
        run.S = run.solv.S
     elif (matrix and multInts):
        run.P = matrix['P']
        
     if not multInts:
        ### calculate multipole integrals 
        run.D,run.Q,run.O,run.H = getM(run.bfs) 
        run.clock.actualize('multipole integrals evaluation')
     else: 
        run.S,run.D,run.Q,run.O,run.H = multInts 

     run.start = True
   else:
     print("Copying the basic variables from singleton object")


