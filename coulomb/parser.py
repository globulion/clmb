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
               exchange=False,
               units='Angstroms'):
      ### --- defaults
      self.units = units

      # ESP options
      self.pot     = pot
      self.svd     = svd
      self.Print   = Print
      self.stat    = stat
      self.pad     = pad
      self.mpot    = mpot 

      # CAMM settings
      self.multints_set = list()
      self.bond_set     = list()

      # Transition density matrix options
      self.transition_set = list()
      self.state_set      = list()
      self.dmat_set       = list()

      # K-integral correction for EEL(1)EDS
      self.exchange = exchange   # ! exchange integrals not working properly!!!
      
      ### --- parse data form input
      self.text = open(file).readlines()
      # derive method, basis set, multipole action
      self.method()
      # withdraw molecule(s) data: name, coords, dmat, smat, multipole integrals
      self.molecules()

  def method(self):
      """Parses method, basis of computation (common 
for all the molecules in the input file!) and 
tasks to be performed"""
      
      first = self.text[0]
      fp = re.compile('# *?([a-zA-Z0-9]*)\/([a-zA-Z0-9]*-?[a-zA-Z0-9|*|+|\(|\)|,]*)\s(.*)')
      fpr= fp.search(first).groups()
      self.method = fpr[0]
      self.basis  = fpr[1]
      self.tasks  = fpr[2].split()
      #self.state  = int(re.search('(?<=state=)[0-9]*', first).group(0))

  def molecules(self):
      """Withdraws molecular specifications:
multiplicity, charge, coordinates, density matrix (optional) 
and multipole integrals (optional)."""
         
      nat = {'H': 1, 'He': 2, 'Li': 3, 'C': 6, 'N': 7, 'O': 8, 'F': 9, 'Na':11, 'S':16, 'Cl':17, 'Mg':12}
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
      for ia, a in enumerate(m):
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

          # read electronic state
          state = 1
          for line in a:
              if line.startswith('STATE='):
                 state = int(line.split('=')[-1])
                 break
          self.state_set.append( state )

          # determine if transition
          trans = False
          for line in a:
              if line.startswith('TRANS='):
                 trans = bool(line.split('=')[-1])
                 break
          self.transition_set.append( trans )

          # read molecular density matrix
          for line in a:
              if line.startswith('DMATRIX='):
                 dmat_file = line.split('=')[1]
                 
                 # read the type of density (SCF, MP2 or CC)
                 for line1 in a:
                     dtype = None
                     if 'DTYPE=' in line1: 
                         dtype = line1.split('=')[-1]
                         break

                 if self.transition_set[ia]:
                    dmat = read_transition_dmatrix(                                             
                               matrix_querry='Alpha transition density to state',
                               filetype='gaussian',file=dmat_file,
                               state=self.state_set[ia])
                                                                                                
                    dmat+= read_transition_dmatrix(
                               matrix_querry='Beta transition density to state',
                               filetype='gaussian',file=dmat_file,
                               state=self.state_set[ia])
                    dmat=sqrt(2.0)*0.5*dmat # due to normalization of coeffs in gaussian to 1/2
                 else:
                    if dtype is not None:
                       dmat = ParseDmatFromFchk(dmat_file,type=dtype)
                    else:
                       print(" No density type specified! Add DTYPE=(SCF, MP2 or CC) keyword.\n Quitting...\n")
                       exit()
                 self.dmat_set.append(dmat)
                 break

          # read connectivity 
          bonds = False
          for line in a:
              if line.startswith('BONDS='):
                 bonds = True
                 self.bond_set.append( [ tuple(array(x.split('-'),int)-1) for x in line.split('=')[-1].split(':') ] )
                 break
          if not bonds: self.bond_set.append( None )

          # read molecular multipole integrals
          for line in a:
              if line.startswith('MULTINTS='):
                 raise NotImplementedError("Reading multipole integrals in not implemented yet! All integrals are as for now evaluated by PyQuante-Mod.")
                 self.multints_set.append([])
                 break
                
      # if not read, please calculate them using PyQuante routines
      if self.multints_set == []:
         pass#for molecule in M:

                
      self.M = M
