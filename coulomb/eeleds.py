# ----------------------------------------------------------- #
#         FIRST-ORDER ELECTROSTATIC INTERACTION ENERGY        #
#            FROM HYBRID VARIATIONAL-PERTURBATIONAL           #
#           INTERACTION ENERGY DECOMPOSITION SCHEME           #
# ----------------------------------------------------------- #

#---------------------------------
from numpy import *
#---------------------------------
import re, os, sys, getopt, random
from time import time
from sys  import argv, exit, stdout
#------------------------------------------------------------------------
from PyQuante.Molecule import Molecule
from PyQuante          import SCF
from PyQuante.Ints     import getbasis, getS, getM, getJ, getK, get2ints, sortints
#from PyQuante.Ints     import sortints
#from PyQuante.cints    import ijkl2intindexBTF as intindexBTF
from PyQuante.CGBF     import coulomb
#-------------------------
from libbbg.units import *
#-------------------------
# memory profiler analyzes memory usage of decorated function
# @memprof(plot = True, threshold = 1024)
#
#from memprof import * 
#-------------------------
# Debug PyQuante routines
#
#import logging
#logging.basicConfig(format="%(message)s",level=logging.DEBUG)
#-------------------------

__all__=['EELEDS']

class EELEDS:
   """contains useful procedures for computing Eint,el,10 from EDS"""
   
   def __init__(self,molecule1,molecule2,basis,method,
                transition=False,exchange=False,
                matrixa=None,matrixb=None):
       self.exchange = exchange
       self.transition = transition
       # form monomers in DCBS
       dimer = self.makeDimer(molecule1,molecule2)
       mon1 = self.makeMonomerDcbs(molecule1,molecule2,2) 
       mon2 = self.makeMonomerDcbs(molecule1,molecule2,1)
       DCBS = getbasis(dimer,basis)
       #bfs1 = getbasis(molecule1,basis)
       #bfs2 = getbasis(molecule2,basis)

       # Minimal sanity check
       y = input("Basis size: %d basis functions.\nTwo-electron integrals will need roughly %4.1f mb.\nWant to continue? (1/0)\n" % (len(DCBS),(len(DCBS)**4)*8/1048576.))
       if not y: exit()
       if self.exchange and not self.transition:
           print "Exchange term is implemented only for transition densities"
           self.exchange = False  
       if self.transition and 'NoneType' in str(type(matrixa)):
           print "Transition densities were not provided in the input file. Error termination."
           exit()
       if not self.transition and not 'NoneType' in str(type(matrixa)):
           print "Transition densities are to be read but the trans keyword is not set in input. Error termination."
           exit()

       if 'NoneType' in str(type(matrixa)):
           # calculate density matrices of monomers in DCBS
           run1 = SCF(molecule1,bfs=DCBS,method=method)
           run1.iterate()
           Pa = 2*run1.dmat
           del run1
           run2 = SCF(molecule2,bfs=DCBS,method=method)
           run2.iterate()
           Pb = 2*run2.dmat
           # save electron repulsion integrals in DCBS
           ERI = run2.ERI
           del run2
       else:
           # read density matrices of monomers in DCBS
           Pa = matrixa
           Pb = matrixb
           # calculate electron repulsion integrals in DCBS
           ERI = get2ints(DCBS)

       # calculate interaction energy!
       self.CalcEint(ERI,DCBS,mon1,mon2,Pa,Pb)

   def makeDimer(self,molecule1,molecule2):
       """makes dimer from molecules 1 and 2"""
       
       position_list = []
       mol_list = [molecule1,molecule2]
       for i in range(2):
           for j in range(len(mol_list[i].atoms)):
               position_list.append(( mol_list[i].atoms[j].atno,
                                      mol_list[i].atoms[j].pos() ))
       mol = Molecule('dimer',position_list,units='Bohr')
       return mol

   def makeMonomerDcbs(self,mol1,mol2,which):
       """sets ghosts bf functions on choosen monomer (which)"""
       
       dimer = self.makeDimer(mol1,mol2)
       N_max = len(dimer.atoms)
       n1 = len(mol1.atoms)
       n2 = len(mol2.atoms)
       if which==1:
          dimer.name = mol2.name
          for i in range(n1):
              dimer.atoms[i].atno = 0
       elif which==2:
          dimer.name = mol1.name
          for i in range(n1,n1+n2):
              dimer.atoms[i].atno = 0
       return dimer

   def getJa(self,bfs,D):
        "Form the Coulomb operator corresponding to a density matrix D"
        nbf = D.shape[0]
        D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
        J = zeros((nbf,nbf),'d')

        for i in xrange(nbf):
            for j in xrange(i+1):

                temp = zeros(nbf*nbf,'d')
                ij = i*(i+1)/2+j
                for k in xrange(nbf):
                    for l in xrange(k+1):
                        kl = k*nbf+l
                        lk = l*nbf+k
                        temp[kl] = coulomb(bfs[i],bfs[j], bfs[k],bfs[l])
                        temp[lk] = temp[kl]

                J[i,j] = dot(temp,D1d)
                J[j,i] = J[i,j]
        return J

   def getJb(self,bfs,D):
        "Form the Coulomb operator corresponding to a density matrix D"
        nbf = D.shape[0]
        D1d = reshape(D,(nbf*nbf,)) #1D version of Dens
        J = zeros((nbf,nbf),'d')
        for i in xrange(nbf):
            for j in xrange(i+1):

                temp = zeros(nbf*nbf,'d')
                kl = 0
                ij = i*(i+1)/2+j
                for k in xrange(nbf):
                    for l in xrange(nbf):
                        temp[kl] = coulomb(bfs[i],bfs[j], bfs[k],bfs[l])
                        kl += 1

                J[i,j] = dot(temp,D1d)
                J[j,i] = J[i,j]
        return J

   #@memprof(plot = True, threshold = 1024)
   def CalcEint(self,ERI,DCBS,mon1,mon2,Pa,Pb):
       """calculates E_el EDS interaction energy"""
       
       Eint = 0
       #K1 = len(MCBS1)
       #Pa1 = self.Block(Pa,0,K1)
       #Pa2 = self.Block(Pa,K1+1,len(Pa))
       #Pb1 = self.Block(Pb,0,K1)
       #Pb2 = self.Block(Pb,K1+1,len(Pa))
       #D1 = self.CalcIntMat(P1b,len(bfs1),len(bfs2),IntsBTF12,1)
       #D2 = self.CalcIntMat(P1a,len(bfs1),len(bfs2),IntsBTF12,2)

       if not self.transition:
          # nuclear attraction terms
          V1 = self.CalcNuclAtt(DCBS,mon2)
          V2 = self.CalcNuclAtt(DCBS,mon1)
          v1b = trace(dot(Pa,V1)) 
          v2a = trace(dot(Pb,V2)) 
          # nuclear repulsion term
          vab = self.CalcNucNuc(mon1,mon2)

       # electron repulsion integrals in-core
       G = getJ(ERI,Pb) 
       if self.exchange:
          X = -0.5*getK(ERI,Pb)
       del ERI
       # direct calculation of electron repulsion integrals
       # getJa should be faster for larger molecules
       #
       #G = self.getJa(DCBS,Pb) 
       #G = self.getJa(DCBS,Pb) 

       # electron repulsion term
       v12 = trace(dot(Pa,G))

       #print "vab = ",vab
       #print "v1b = ",v1b
       #print "v2a = ",v2a
       #print "v12 = ",v12

       if not self.transition:
          self.eint_coul = (vab+v1b+v2a+v12) * UNITS.HartreeToCmRec
       else:
          self.eint_coul = v12 * UNITS.HartreeToCmRec
       if self.exchange: 
          self.eint_exch=trace(dot(Pa,X)) * UNITS.HartreeToCmRec

   def CalcNucNuc(self,MOL1,MOL2):
       """calculates interaction energy between two sets of nuclei 
          in DCBS molecule objects"""
          
       Eint = 0
       for i in MOL1.atoms:
           for j in MOL2.atoms:
               if i.atid!=j.atid:
                  R = sqrt(sum(( array(i.pos())-array(j.pos()) )**2))
                  Eint+= i.atno*j.atno / R
       return Eint

   def CalcNuclAtt(self,bfs_A,mol_B):
       """calculates attraction energy matrix between electrons 
          from molecule A and nuclei from molecule B"""
          
       K = len(bfs_A)
       V = zeros((K,K),dtype=float64)
       for i in xrange(K):
           for j in xrange(K):
               for atom in mol_B.atoms:
                   V[i,j]+= atom.atno* bfs_A[i].nuclear(bfs_A[j],atom.pos())
       return V

   def Block(self,matrix,K1,K2):
        """makes a block from square matrix between K1 and K2 column"""
        
        Ktot = len(matrix)
        m = zeros((Ktot,Ktot),dtype=float64)
        m[K1:K2,K1:K2] = matrix[K1:K2,K1:K2]
        return m

   def __repr__(self):
       """prints information about interaction energy result"""
       log = "\n"
       log+= " ------------------------------------\n"
       log+= " INTERACTION ENERGY COMPONENTS [CM-1]\n"
       log+= " ------------------------------------\n"
       log+=    " COULOMB    : %12.6f\n" %  self.eint_coul
       if self.exchange:
          log+= " EXCHANGE   : %12.6f\n" %  self.eint_exch
          log+= " TOTAL      : %12.6f\n" % (self.eint_coul+self.eint_exch)
       else:
          log+= " ECHANGE REPULSION TERM NOT ESTIMATED\n"
       log+= " ------------------------------------\n"
       
       return str(log)
