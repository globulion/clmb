# --------------------------------------------- #
#         ELECTROSTATIC POTENTIAL FITTING       #
# --------------------------------------------- #

from run            import *
from multip         import *
from elpot          import *
from libbbg.units   import *

__all__ = ['ESP']

class ESP(RUN):
   """contains usefull procedures for fitting charges 
based on density matrix and potential. Potential can 
be calculated from CAMMs or directly from the wave 
function ('WFN'). The default is CAMM. Can perform 
SVD decomposition of distance matrix in order to 
verify the quality of charge fitting. All values in AU."""

   def __init__(self, molecule, basis, method, mpot=1000, pot='CAMM', pad=10, 
                      SVD=False, matrix=None, multInts=None, cnt=0.0100, 
                      stat=False, Print=True, transition=False):
       RUN.__init__(self, molecule, basis, method, matrix, multInts, hexadecapoles=False)
       # transition 
       self.transition = transition
       # number of points for potential to be computed (total)
       self.mpot = mpot * len(molecule.atoms)
       self.pot = pot
       # box padding
       self.pad = pad
       # set vdW radii for nuclei atomic number (in au)
       self.vdW = {1:1.5000,
                   6:3.0000, 
                   7:2.4000,
                   8:2.6000,
                   9:2.3000,}
       # condition number treshold
       self.condition_number_treshold = cnt
       self.stat=False
       self.svd =False
       # --- main action ---
       # calculate potential in some points
       self.RandomPoints()
       self.CalcPot()
       # svd analysis
       if SVD: 
          self.svd = True
          self.rijmatrix()
       # final fitting
       self.charges = self.Final_fit(self.lsmatrix(),self.bvec())[:-1]
       # print the results
       if Print: 
          self.__print__()
       # check the quality of fitting
       if stat:
          self.stat=True
          self.CalcPot_from_charges()
          self.statistics()

   def Final_fit(self,A,B):
       """performs final fitting of charges"""
       return dot(linalg.inv(A),B)

   def CalcPot(self):
       """calculates potential in given points""" 
       
       # derive potential from CAMMs
       if   self.pot.lower()=='camm':
            camm = MULTIP(self.molecule,self.basis,self.method,
                          matrix=self.matrix,multInts=self.multInts,
                          transition=self.transition)
            camm.camms()
            self.clock.actualize('calculate CAMMs for potential estimation')
            camm.makeTracelessCAMMs()
            self.clock.actualize('convert CAMMs to traceless tensors')
            for i in range(self.mpot):
                self.VArray[i,3] = Vr_camms(camm,self.VArray[i,:3]) 
            self.clock.actualize('potential from CAMMs')
       elif self.pot.lower()=='wfn':
            for i in range(self.mpot):
                self.VArray[i,3] = Vr_wfn_1(self.molecule,self.bfs,self.P,self.VArray[i,:3])
            self.clock.actualize('potential from wave function')
       else:print "\nnot yet written. Quitting COULOMB.py ...\n"; exit()

   def CalcPot_from_charges(self):
       """calculates the potential from fitted charges"""
       
       VFArray = zeros((self.mpot,4),dtype=float64)
       VFArray[:,:3] = self.VArray[:,:3]
       for i in range(self.mpot):
           for j in range(len(self.charges)):
               VFArray[i,3] += self.charges[j]/sqrt( sum( (VFArray[i,:3]-self.RArray[j])**2 ) )
       self.VFArray = VFArray
       self.clock.actualize('potential from fitted charges')

   def __print__(self):
       """an output printout routine. Prints: 
   --- real order of distance matrix (if specified)
   --- fitted charges
   --- statistical assessment of fitting quality (if specified)"""
   
       log = '\n --- ESP procedure ---\n'
       if self.svd: 
          log+= "\n     THE ORDER OF DISTANCE MATRIX IS --- %d ---"%self.order
          log+= "\n     for condition number treshold: %.6f\n"%self.condition_number_treshold
          if self.order<self.N_at: 
             log+= "\n     LINEAR DEPENDENCE IN DISTANCE MATRIX DETECTED! "
             log+= "\n     --- %d poitns were discarded from the analysis"%(self.N_at-self.order)
          else:
             log+= "\n     NO LINEAR DEPENDENCE IN DISTANCE MATRIX."

       log+= '\n\n\n\n'
       log+= " FITTED CHARGES\n"
       log+= " --------------\n"
       log+= " %s %s\n"                         % ('ATOM'.rjust(12), 'q'.rjust(12))
       for i in range(self.N_at):
           log+= " %10d %12.8f\n"               % ( i, self.charges[i] )  
       log+= " \n"   
       print log

   # ---- UTILITIES
   def RandomPoints(self):
       """compute points for potential in a random manner. In au."""
       
       # find cube of molecule 
       xmin,xmax = min(self.RArray[:,0]),max(self.RArray[:,0])
       ymin,ymax = min(self.RArray[:,1]),max(self.RArray[:,1])
       zmin,zmax = min(self.RArray[:,2]),max(self.RArray[:,2])
       MIN = array([xmin,ymin,zmin]) - self.pad
       MAX = array([xmax,ymax,zmax]) + self.pad
       BOX = MAX-MIN
       # random points of potential
       VArray = zeros((self.mpot,4),dtype=float64)
       i=0
       while i<self.mpot:
           r=0
           X = MIN[0] + random.random() * BOX[0]
           Y = MIN[1] + random.random() * BOX[1]
           Z = MIN[2] + random.random() * BOX[2]
           p = array([X,Y,Z,0])
           # check if the point lies outside the vdW space
           for j in range(self.N_at):
               r += self.vdW_sphere(self.vdW[self.molecule.atoms[j].atno],p[:3],self.RArray[j])
           if not r: 
              VArray[i,:] = p
              i+=1
           else: continue
       self.clock.actualize('generation of %d random points' % self.mpot )
       self.VArray = VArray

   def rijmatrix(self):
       """calculates distance matrix, A_ij = 1./r_ij"""
       
       A = zeros((self.mpot,self.N_at),dtype=float64)
       for i in range(self.mpot):
           for j in range(self.N_at):
               A[i,j] = 1./sqrt( sum( (self.VArray[i,:3]-self.RArray[j])**2 ) )
       self.clock.actualize('1./rij matrix calculations')
       # check the order of the matrix
       U,S,Vt = svd(A)
       del U
       del Vt
       self.clock.actualize('SVD decomposition of distance matrix')
       Smax = max(S)
       order = self.N_at
       for i in range(len(S)):
           if S[i]/Smax < self.condition_number_treshold:
              order -=1
       del S
       self.order = order
       pass
       #return A,order

   def lsmatrix(self):
       """computes least-square matrix, A_jk=\sum_i (r_ij*r_ik)^-1"""
       
       A = zeros((self.N_at+1,self.N_at+1),dtype=float64)
       for j in range(self.N_at):
           for k in range(j+1):
               for i in range(self.mpot):
                   rij = sqrt( sum( (self.VArray[i,:3]-self.RArray[j])**2 ) )
                   rik = sqrt( sum( (self.VArray[i,:3]-self.RArray[k])**2 ) )
                   A[j,k] += 1./(rij*rik)
               A[k,j] = A[j,k]
       A[-1,-1] = 1
       self.clock.actualize('computation of least-square matrix')
       return A

   def bvec(self):
       """computes B vector, B_k = \sum_i V_i/r_ik"""
       
       B = zeros(self.N_at+1,dtype=float64)
       for k in range(self.N_at):
           for i in range(self.mpot):
               B[k] += self.VArray[i,3]/sqrt( sum( (self.VArray[i,:3]-self.RArray[k])**2 ) ) 
       self.clock.actualize('B-matrix computation')
       return B

   def vdW_sphere(self,vdW_radius,R,Ri):
       """checks if point lies within atom's vdW sphere."""
       """if '0' - the point lies outside atom's vdW sphere."""
       
       r2 = sum((R-Ri)**2,axis=0)
       if r2>vdW_radius**2: return 0
       else:                return 2

   def statistics(self):
       """makes statistical assessment of the quality of charges fitting
          basing on the potential calculated from fitted charges and reference potential"""
          
       V     = self.VFArray[:,-1]
       V_ref = self.VArray[:,-1]
       import matplotlib.pyplot as plt
       # relative error in per cent unit
       rel_err = (V-V_ref)/V_ref*100
       plt.hist(rel_err,bins=500,range=(-100,100),histtype='stepfilled')
       plt.show()
       log = '\n --- STATISTICAL EVALUATION OF THE FITTING QUALITY ---\n'
       log+= '\nnot yet written\n'
       print log 
