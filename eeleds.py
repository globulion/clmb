# ----------------------------------------------------------- #
#         FIRST-ORDER ELECTROSTATIC INTERACTION ENERGY        #
#            FROM HYBRID VARIATIONAL-PERTURBATIONAL           #
#           INTERACTION ENERGY DECOMPOSITION SCHEME           #
# ----------------------------------------------------------- #

from coulomb_head import *
from units        import *

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
       print "Basis size: %d basis functions. Want to continue?" % len(DCBS)
       y = input("1/0")
       if not y: exit()
       #bfs1 = getbasis(molecule1,basis)
       #bfs2 = getbasis(molecule2,basis)
       # calculate density matrices of monomers in DCBS
       run1 = SCF(molecule1,bfs=DCBS,method=method)
       run1.iterate()
       if 'NoneType' in str(type(matrixa)):
           Pa = 2*run1.dmat
           run2 = SCF(molecule2,bfs=DCBS,method=method)
           run2.iterate()
           Pb = 2*run2.dmat
       else:
           Pa = matrixa
           Pb = matrixb
       # 2-electron integrals 
       ERI = run1.ERI
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


   def CalcEint(self,ERI,DCBS,mon1,mon2,Pa,Pb):
       """calculates E_el EDS interaction energy"""
       
       Eint = 0
       #K1 = len(MCBS1)
       #Pa1 = self.Block(Pa,0,K1)
       #Pa2 = self.Block(Pa,K1+1,len(Pa))
       #Pb1 = self.Block(Pb,0,K1)
       #Pb2 = self.Block(Pb,K1+1,len(Pa))

       if not self.transition:
          V1 = self.CalcNuclAtt(DCBS,mon2)
          V2 = self.CalcNuclAtt(DCBS,mon1)
       #D1 = self.CalcIntMat(P1b,len(bfs1),len(bfs2),IntsBTF12,1)
       #D2 = self.CalcIntMat(P1a,len(bfs1),len(bfs2),IntsBTF12,2)
       G = getJ(ERI,Pb) 
       if self.exchange: X = -0.5*getK(ERI,Pb)
       if not self.transition:
          vab = self.CalcNucNuc(mon1,mon2)
          v1b = trace(dot(Pa,V1)) 
          v2a = trace(dot(Pb,V2)) 
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