# ------------------------------------------- #
#         MULTIPOLE DISTRIBUTION MODULE       #
# ------------------------------------------- #

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
from util import read_transition_dmatrix
# ---------------------

__all__=['PARSER']

class PARSER:
  """withdraws data from input file"""
  
  def __init__(self,file,
               pot='CAMM',
               svd=False,
               Print=True,
               stat=False,
               pad=10,
               mpot=1000,
               state=1,
               transition=False,
               exchange=False,
               units='Angstroms',
               dmat_file_type='gaussian'):
      ### --- defaults
      self.units = units
      # density matrix file type
      self.dmat_file_type = dmat_file_type

      # ESP options
      self.pot = pot
      self.svd = svd
      self.Print = Print
      self.stat = stat
      self.pad = pad
      self.mpot = mpot 
      
      # Transition density matrix options
      self.transition = transition
      self.state = state
      # K-integral correction for EEL(1)EDS
      self.exchange = exchange
      
      ### --- parse data form input
      self.text = open(file).readlines()
      # derive method, basis set, multipole action
      self.method()
      # withdraw molecule(s) data: name, coords, dmat, smat, multipole integrals
      self.molecules()

  def method(self):
      """parses method, basis of computation (common 
      for all the molecules in the input file!) and 
      tasks to be performed"""
      
      first = self.text[0]
      fp = re.compile('# *?([a-zA-Z0-9]*)\/([a-zA-Z0-9]*-?[a-zA-Z0-9/**/(*)*]*)\s(.*)')
      fpr= fp.search(first).groups()
      self.method = fpr[0]
      self.basis  = fpr[1]
      self.tasks  = fpr[2].split()
      self.state  = int(re.search('(?<=state=)[0-9]*', first).group(0))

  def molecules(self):
      """withdraws molecular specifications:
         multiplicity, charge, coordinates, density matrix (optional) 
         and multipole integrals (optional)."""
         
      self.dmat_set = []
      self.multints_set = [] 
      nat = {'H': 1, 'He': 2, 'Li': 3, 'C': 6, 'N': 7, 'O': 8, 'F': 9}
      # search atomic coorinates
      cp = re.compile('([a-zA-Z]*) *(-?\d*\.\d*) *(-?\d*\.\d*) *(-?\d*\.\d*).*')
      # search multiplicity and charge
      mc = re.compile('^(-?\d) *(\d).*')
      Y=[]
      for i in range(len(self.text)):
          if self.text[i].startswith('&MOL'): Y.append(i)

      m = []
      M = []
      i=1
      for mol in Y:
          y=self.text[mol+i]
          r=self.text[mol+i]
          while not y.startswith('&'):
                i+=1
                y =self.text[mol+i]
                r+=y
          m.append(r.split('\n'))
          i=1
      # read molecular data
      for a in m:
          name = a[0]
          coords= []
          # multiplicity and charge
          for line in a:
              result = mc.search(line)
              if result: 
                 t = result.groups()
                 charge = int(t[0])
                 multiplicity = int(t[1])
          # coords
          for line in a:
              result = cp.search(line)
              if result:
                 t = result.groups()
                 coord = ( nat[t[0]], ( float(t[1]), float(t[2]), float(t[3]) ) )
                 coords.append(coord)
          mol1 = Molecule('%s' % name,coords,multiplicity=multiplicity,charge=charge,units=self.units)
          M.append(mol1)
          # try to read molecular density matrix and overlap matrix 
          for line in a:
              if line.startswith('DMATRIX='):
                 dmat = read_transition_dmatrix(
                            matrix_querry='Alpha transition density to state',
                            filetype=self.dmat_file_type,file=line.split('=')[-1],
                            state=self.state)
                
                 dmat+= read_transition_dmatrix(
                            matrix_querry='Beta transition density to state',
                            filetype=self.dmat_file_type,file=line.split('=')[-1],
                            state=self.state)
                 dmat=sqrt(2.0)*0.5*dmat # due to normalization of coeffs in gaussian to 1/2

                 self.dmat_set.append(dmat)
                 break
          # read molecular multipole integrals
          for line in a:
              if line.startswith('MULTINTS='):
                 self.multints_set.append([])
                 break
                
      # if not read, please calculate them using PyQuante routines
      if self.multints_set == []:
         pass#for molecule in M:

                
      self.M = M
