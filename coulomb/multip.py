# ------------------------------------------- #
#         MULTIPOLE DISTRIBUTION MODULE       #
# ------------------------------------------- #

#from libbbg.dma import DMA
#from run        import *
#from libbbg.utilities2 import numpy.array_outer_product,    \
#                              numpy.array_outer_product_1_2,\
#                              numpy.array_outer_product_2_1
import numpy, run, libbbg.utilities, libbbg.dma

__all__=['MULTIP']

class MULTIP(run.RUN):
    """\                                                          
==================================================================
Ccontains useful procedures for computing:                        
  - Molecular Multipole Moments (MMMs)                            
  - Cumulative Atomic Multipole Moments (CAMMs)                   
  - Cumulative Atomic and Bond Multipole Moments (CABMMs)         
------------------------------------------------------------------
USAGE:                                                            
  <object>=MULTIP(molecule,basis,method,matrix=None,multInts=None,
                  transition=False,bonds=[],hexadecapoles=False)                      
------------------------------------------------------------------
1) calculating the moments:                                       
  <object>.mmms()                                                 
  <object>.camms()                                                
2) returns the list of DMA objects created in <object>            
  <object>.get()                                                  
3) old printouts                                                  
  <object>.__printMMMs__()                                        
  <object>.__printCAMMs__()                                       
4) new printout                                                   
  print <object>                                                  
==================================================================
"""
       
    def __init__(self, molecule, basis, method,
                       matrix=None,multInts=None,transition=False,
                       bonds=None,vec=None,hexadecapoles=False):
        run.RUN.__init__(self, molecule, basis, method, matrix, multInts, hexadecapoles)
        # LIST1 - the list of atoms in the order of basis functions used,
        # e.g.: for h2o molecule with atoms: 8,1,1 and STO-3G basis
        # LIST1 = 0    0    0    0    0    1    2
        #         |    |    |    |    |    |    |
        #        1s   2s   2px  2py  2pz  1s   1s
        self.LIST1 = self.bfs.LIST1
        self.transition = transition
        self.name=molecule.name
        self.bonds = bonds
        self.vec = vec
        # bonds - the list of tuples of atomic pairs engaged in bonds:
        # bonds = [ (n11,n12), (n21,n22), ... , (ni1,ni2) ]           
        # where ni1 > ni2 specifies list of i bonds, n are indices    
        # in Python convention (N-1)                                  
        self.pos, self.origin = self.__make_pos()
        self.__dma_bin = []
        self.operation = 'None'
        self.has_hexadecapoles = hexadecapoles

    def get(self):
        """returns the DMA objects from all the runs"""
        return self.__dma_bin
    
    def mmms(self):
        """calculate molecular multipole moments"""
        
        self.DipoleMoment    = self.MU()
        self.QuarupoleMoment = self.QUAD()
        self.OctupoleMoment  = self.OCT()
        if self.has_hexadecapoles: 
           self.HexadecapoleMoment = self.HEX()
        self.operation = 'MMM'
        self.clock.actualize("calculation of MMMs")
        # save the DMA object
        #moments = [0,self.DipoleMoment,self.QuarupoleMoment,self.OctupoleMoment]
        #self.__make_dma(moments=moments,
        #                change_origins=False)  # napraw to bo to nie dziala!
        return
        
    def get_centroids(self):
        """calculate LMO centroids"""
        nmos = len(self.vec)
        centroids = numpy.zeros((nmos,3),dtype=numpy.float64)
        for LMO in xrange(nmos):
            centroids[LMO] = numpy.tensordot( self.vec[LMO], numpy.tensordot(self.vec[LMO],
                                                                 self.D,(0,1)), (0,1))
        return centroids
    
    def camms(self):
        """evaluates the distributed moments"""
        if self.bonds is not None: self.__cabmms()
        elif self.vec is not None: self.__lmtp()
        else:                      self.__camms()
        return


    def __lmtp(self):
        """calculate LMTP distribution"""
        self.operation = 'LMTP'
        nmos = len(self.vec)
        
        Mon      = []
        Dip      = []
        Quad     = []
        Oct      = []
        Hex      = []
        
        # nuclear moments
        for atom in self.molecule.atoms:
            if self.transition: qA  = 0
            else:               qA  = atom.atno
            qA = atom.atno
            R  = numpy.array(atom.pos())
            RR = numpy.outer(R,R)
            RRR= numpy.outer(R,numpy.outer(R,R)).reshape(3,3,3)
            MA = qA * R
            QA = qA * RR
            OA = qA * RRR
            #
            Mon .append(qA)
            Dip .append(MA)
            Quad.append(QA)
            Oct .append(OA)
            if self.has_hexadecapoles:
               HA = qA * numpy.outer(R,RRR).reshape(3,3,3,3)
               Hex.append(HA)
            
        # electronic moments
        for LMO in xrange(nmos):
            A = self.vec[LMO]
            chr = - 2* numpy.tensordot( A, numpy.tensordot(A, self.S,(0,0)), (0,0))
            dip = - 2* numpy.tensordot( A, numpy.tensordot(A, self.D,(0,1)), (0,1))
            qdr = - 2* numpy.tensordot( A, numpy.tensordot(A, self.Q,(0,2)), (0,2))
            oct = - 2* numpy.tensordot( A, numpy.tensordot(A, self.O,(0,3)), (0,3))
            #
            Mon .append(chr)
            Dip .append(dip)
            Quad.append(qdr)
            Oct .append(oct)
            if self.has_hexadecapoles:
               hex = - 2* numpy.tensordot( A, numpy.tensordot(A, self.H,(0,4)), (0,4))
               Hex.append(hex)


        # save LMTP
        self.Mon   = Mon
        self.Dip   = Dip
        self.Quad  = Quad
        self.Oct   = Oct
        if self.has_hexadecapoles:
           self.Hex = Hex
        
        # clock measure
        self.clock.actualize('computing LMTP')
        
        # create the DMA object for calculated property
        self.__make_dma(change_origins=True)
        return
        
    def __cabmms(self):
        """calculates C(A+B)MMs"""
        self.operation = 'CABMM'
        
        Mon      = []
        Dip      = []
        Quad     = []
        Oct      = []
        Hex      = []
        
        n = 0 # atom number
        ### [1] atomic moments
        for atom in self.molecule.atoms:
            if self.transition: qA  = 0
            else:               qA  = atom.atno
            R  = numpy.array(atom.pos())
            RR = numpy.outer(R,R)
            RRR= numpy.outer(R,numpy.outer(R,R)).reshape(3,3,3)
            MA = qA * R
            QA = qA * RR
            OA = qA * RRR
            if self.has_hexadecapoles:
               HA = qA * numpy.outer(R,RRR).reshape(3,3,3,3)
          
            for I in xrange(self.K):
                i = self.LIST1[I]
                for J in xrange(self.K):
                    j =  self.LIST1[J]
                    if ((i==n==j) or ((i==n and j!=i) and (not self.__IJ_are_bonded(i,j))) ):
                       qA -= self.P[I,J] * self.S[I,J]
                       MA -= self.P[I,J] * self.D[:,I,J]
                       QA -= self.P[I,J] * self.Q[:,:,I,J]
                       OA -= self.P[I,J] * self.O[:,:,:,I,J]
                       if self.has_hexadecapoles:
                          HA -= self.P[I,J] * self.H[:,:,:,:,I,J]
                    #elif ((i==n and j!=i) and (not self.__IJ_are_bonded(i,j))):
                    #   qA -= self.P[I,J] * self.S[I,J]
                    #   MA -= self.P[I,J] * self.D[:,I,J]
                    #   QA -= self.P[I,J] * self.Q[:,:,I,J]
                    #   OA -= self.P[I,J] * self.O[:,:,:,I,J]

            Mon .append(qA)
            Dip .append(MA)
            Quad.append(QA)
            Oct .append(OA)
            if self.has_hexadecapoles:
               Hex .append(HA)
            n+=1
            
        ### [2] bond moments
        # code for python-2.7
        qB  = {bond:0                            for bond in self.bonds}
        MB  = {bond:numpy.zeros( 3     ,dtype=numpy.float64) for bond in self.bonds}
        QB  = {bond:numpy.zeros((3,3  ),dtype=numpy.float64) for bond in self.bonds}
        OB  = {bond:numpy.zeros((3,3,3),dtype=numpy.float64) for bond in self.bonds}
        if self.has_hexadecapoles:
           HB  = {bond:numpy.zeros((3,3,3,3),dtype=numpy.float64) for bond in self.bonds}
        # code for python-2.6
        #qB = {}
        #MB = {}
        #QB = {}
        #OB = {}
        #HB = {}
        #for bond in self.bonds:
        #    qB.update({bond:0})
        #    MB.update({bond:numpy.zeros( 3     ,dtype=numpy.float64)})
        #    QB.update({bond:numpy.zeros((3,3  ),dtype=numpy.float64)})
        #    OB.update({bond:numpy.zeros((3,3,3),dtype=numpy.float64)})
        #    if self.has_hexadecapoles:
        #       HB.update({bond:numpy.zeros((3,3,3,3),dtype=numpy.float64)})
        #
        for bond in self.bonds:
            for I in xrange(self.K):
                for J in xrange(self.K):
                    i = self.LIST1[I]
                    j = self.LIST1[J]
                    if (i,j) == bond:
                        qB[bond] -= 2*self.P[I,J] * self.S[I,J]
                        MB[bond] -= 2*self.P[I,J] * self.D[:,I,J]
                        QB[bond] -= 2*self.P[I,J] * self.Q[:,:,I,J]
                        OB[bond] -= 2*self.P[I,J] * self.O[:,:,:,I,J]
                        if self.has_hexadecapoles:
                           HB[bond] -= 2*self.P[I,J] * self.H[:,:,:,:,I,J]
        # append the bond moments
        for bond in self.bonds:
            Mon .append( qB[bond] )
            Dip .append( MB[bond] )
            Quad.append( QB[bond] )
            Oct .append( OB[bond] )
            if self.has_hexadecapoles:
               Hex .append( HB[bond] )
        # save C(A+B)MMS 
        self.Mon   = Mon
        self.Dip   = Dip
        self.Quad  = Quad
        self.Oct   = Oct
        if self.has_hexadecapoles:
           self.Hex = Hex

        # clock measure
        if self.bonds: self.clock.actualize('computing C(A+B)MMs')
        else:          self.clock.actualize('computing CAMMs')
        
        # create the DMA object for calculated property
        self.__make_dma(change_origins=True)
        return
                
    def __camms(self, old_code=False):
        """calculate CAMMs"""
        self.operation = 'CAMM'

        Mon      = []
        Dip      = []
        Quad     = []
        Oct      = []
        Hex      = []
        n = 0    # atom number

        if old_code:
           for atom in self.molecule.atoms:                                                                                                            
               if self.transition: q  = 0
               else:               q  = atom.atno
               M  = numpy.zeros( 3,dtype=numpy.float64)
               Q  = numpy.zeros((3,3),dtype=numpy.float64)
               O  = numpy.zeros((3,3,3),dtype=numpy.float64)
               if self.has_hexadecapoles:
                  H = numpy.zeros((3,3,3,3),dtype=numpy.float64)
               R = numpy.array(atom.pos()) 
               for I in xrange(self.K):
                   if self.LIST1[I]==n:
                      for J in xrange(self.K):
                          q -= self.P[I,J] * self.S[I,J]
                          M += self.P[I,J] * ( self.S[I,J] * R - self.D[:,I,J] )
                          for i in [0,1,2]:
                              for j in [0,1,2]:
                                  Q[i,j] += self.P[I,J] * ( R[i] * self.D[j,I,J] + R[j] * self.D[i,I,J] -R[i] * R[j] * self.S[I,J] - self.Q[i,j,I,J] )
                          #Q += self.P[I,J] * (\
                          #                   outer(R,self.D[:,I,J]) + transpose(outer(R,self.D[:,I,J]),axes=(1,0)) \
                          #                 - outer(R,R) * self.S[I,J] - self.Q[:,:,I,J]                  \
                          #                   )
                          for i in [0,1,2]:
                              for j in [0,1,2]:
                                  for k in [0,1,2]:
                                      O[i,j,k] += self.P[I,J] * ( R[i] * R[j] * R[k] * self.S[I,J] - R[i] * R[j] * self.D[k,I,J] 
                                                              - R[i] * R[k] * self.D[j,I,J] - R[k] * R[j] * self.D[i,I,J] 
                                                              + R[i] * self.Q[j,k,I,J] + R[j] * self.Q[k,i,I,J] + R[k] * self.Q[i,j,I,J] 
                                                              - self.O[i,j,k,I,J] )
                          #R2 = outer(R,R)
                          #R2D= numpy.array_outer_product_2_1(R2,self.D[:,I,J])
                          #RQ = numpy.array_outer_product_1_2(R,self.Q[:,:,I,J])
                          #O += self.P[I,J] * (\
                          #                   numpy.array_outer_product_1_2(R,outer(R,R)) * self.S[I,J] \
                          #                 - R2D - transpose(R2D,axes=(0,2,1)) - transpose(R2D,axes=(2,1,0))\
                          #                 + RQ  + transpose(RQ ,axes=(1,2,0)) + transpose(RQ ,axes=(2,0,1))\
                          #                 - self.O[:,:,:,I,J] 
                          #                   )  
                          if self.has_hexadecapoles:
                             for i in [0,1,2]:
                                 for j in [0,1,2]:
                                     for k in [0,1,2]:
                                         for l in [0,1,2]:
                                             H[i,j,k,l] += self.P[I,J] * (- R[i] * R[j] * R[k] * R[l] * self.S[        I,J] 
                                                                          + R[i] * R[j] * R[k]        * self.D[l,      I,J]
                                                                          + R[i] * R[j] * R[l]        * self.D[k,      I,J]
                                                                          + R[j] * R[k] * R[l]        * self.D[i,      I,J]
                                                                          + R[i] * R[k] * R[l]        * self.D[j,      I,J]
                                                                          - R[i] * R[j]               * self.Q[k,l,    I,J]
                                                                          - R[i] * R[k]               * self.Q[j,l,    I,J]
                                                                          - R[i] * R[l]               * self.Q[j,k,    I,J]
                                                                          - R[j] * R[l]               * self.Q[i,k,    I,J]
                                                                          - R[j] * R[k]               * self.Q[i,l,    I,J]
                                                                          - R[k] * R[l]               * self.Q[i,j,    I,J]
                                                                          + R[i]                      * self.O[j,k,l,  I,J]
                                                                          + R[j]                      * self.O[i,k,l,  I,J]
                                                                          + R[k]                      * self.O[i,j,l,  I,J]
                                                                          + R[l]                      * self.O[i,j,k,  I,J]
                                                                          -                             self.H[i,j,k,l,I,J])
                                                                                                                                                       
               Mon .append(q)
               Dip .append(M)
               Quad.append(Q)
               Oct .append(O)
               if self.has_hexadecapoles:
                  Hex .append(H)
               n+=1
        # newer, much much faster code similar to old CABMM code
        else:  
            for atom in self.molecule.atoms:
                if self.transition: qA  = 0
                else:               qA  = atom.atno
                R  = numpy.array(atom.pos())
                RR = numpy.outer(R,R)
                RRR= numpy.outer(R,numpy.outer(R,R)).reshape(3,3,3)
                MA = qA * R
                QA = qA * RR
                OA = qA * RRR
                if self.has_hexadecapoles:
                   HA = qA * numpy.outer(R,RRR).reshape(3,3,3,3)

                for I in xrange(self.K):
                    i = self.LIST1[I]
                    if i==n:
                       for J in xrange(self.K): 
                           j =  self.LIST1[J]

                           qA -= self.P[I,J] * self.S[I,J]
                           MA -= self.P[I,J] * self.D[:,I,J]
                           QA -= self.P[I,J] * self.Q[:,:,I,J]
                           OA -= self.P[I,J] * self.O[:,:,:,I,J]
                           if self.has_hexadecapoles:
                              HA -= self.P[I,J] * self.H[:,:,:,:,I,J]
                                                                                                  
                Mon .append(qA)
                Dip .append(MA)
                Quad.append(QA)
                Oct .append(OA)
                if self.has_hexadecapoles:
                   Hex .append(HA)
                n+=1
           
        # save CAMMS 
        self.Mon   = Mon
        self.Dip   = Dip
        self.Quad  = Quad
        self.Oct   = Oct
        if self.has_hexadecapoles:
           self.Hex = Hex
 
        # clock measure
        self.clock.actualize('computing CAMMs')

        # create the DMA object for calculated property
        if old_code: self.__make_dma(change_origins=False)
        else:        self.__make_dma(change_origins=True )
        return

    def __IJ_are_bonded(self,I,J):
        """\                                   
determines wheather atoms I and J              
are bound to each other. Returns boolean value 
basing on the self.bonds list of bonds.        
"""
        is_bonded = False
        if I!=J:
           for bond in self.bonds:
               if (I in bond) and (J in bond):
                   is_bonded = True
                   break
        return is_bonded
    
    def __make_pos(self):
        """calculate positions and origins"""
        natoms = len(self.molecule.atoms)
        pos    = numpy.zeros((len(self.molecule.atoms),3),dtype=numpy.float64)
        ### CBAMM
        if self.bonds is not None:
           nfrag  = natoms + len(self.bonds)
           origin = numpy.zeros((nfrag,3),dtype=numpy.float64)
           for i,atom in enumerate(self.molecule.atoms):
               pos[i] = numpy.array(atom.pos())
               origin[i] = numpy.array(atom.pos())
           for i,bond in enumerate(self.bonds):
               A = numpy.array(self.molecule.atoms[bond[0]].pos())
               B = numpy.array(self.molecule.atoms[bond[1]].pos())
               origin[natoms+i] = 0.5*(A+B)
        ### CAMM and MMM
        elif (self.bonds is None) and (self.vec is None):
           nfrag  = natoms
           origin = numpy.zeros((nfrag,3),dtype=numpy.float64)
           for i,atom in enumerate(self.molecule.atoms):
               pos[i] = numpy.array(atom.pos())
               origin[i] = numpy.array(atom.pos())            
        ### LMTP
        elif self.vec is not None:
           self.centroids = self.get_centroids()
           nfrag  = natoms + len(self.centroids)
           origin = numpy.zeros((nfrag,3),dtype=numpy.float64)
           for i,atom in enumerate(self.molecule.atoms):
               pos[i] = numpy.array(atom.pos())
               origin[i] = numpy.array(atom.pos())
           for i,centroid in enumerate(self.centroids):
               origin[natoms+i] = self.centroids[i]
            
        return pos, origin

    def __make_dma(self,moments=None,change_origins=False):
        """creates the DMA object and saves it into self.dma"""
        #if moments is None:
        #   mon=numpy.array(self.Mon),
        #   dip=numpy.array(self.Dip),
        #   quad=numpy.array(self.Quad),
        #   oct=numpy.array(self.Oct)
        #else:
        #   print "FFF"
        #   mon,dip,quad,oct = moments
        #   mon = numpy.array(mon)[newaxis]
        #   dip = numpy.array(dip)[newaxis]
        #   quad= numpy.array(quad)[newaxis]
        #   oct = numpy.array(oct)[newaxis]
        result = libbbg.dma.DMA(nfrag=len(self.origin), hexadecapoles=self.has_hexadecapoles)
        result.set_name("%s --- %s %s/%s" % (self.molecule.name, self.operation,
                                             self.method, self.basis) )
        # in the case of C(A+B)MMs
        if change_origins:
           result.set_structure(pos=self.pos,origin=numpy.zeros((len(self.origin),3),dtype=numpy.float64))
        # in the case of CAMMs and MMMs
        else:
           result.set_structure(pos=self.pos,equal=True)
           
        # accumulate the distributed moments to the DMA object
        ###print numpy.array(self.Quad).shape,numpy.array(quad).shape
        result.DMA[0][:] = numpy.array(self.Mon)
        #
        result.DMA[1][:] = numpy.array(self.Dip)
        #
        result.DMA[2][:,0] = numpy.array(self.Quad)[:,0,0]
        result.DMA[2][:,1] = numpy.array(self.Quad)[:,1,1]
        result.DMA[2][:,2] = numpy.array(self.Quad)[:,2,2]
        result.DMA[2][:,3] = numpy.array(self.Quad)[:,0,1]
        result.DMA[2][:,4] = numpy.array(self.Quad)[:,0,2]
        result.DMA[2][:,5] = numpy.array(self.Quad)[:,1,2]
        #
        result.DMA[3][:,0] = numpy.array(self.Oct)[:,0,0,0]
        result.DMA[3][:,1] = numpy.array(self.Oct)[:,1,1,1]
        result.DMA[3][:,2] = numpy.array(self.Oct)[:,2,2,2]
        result.DMA[3][:,3] = numpy.array(self.Oct)[:,0,0,1]
        result.DMA[3][:,4] = numpy.array(self.Oct)[:,0,0,2]
        result.DMA[3][:,5] = numpy.array(self.Oct)[:,0,1,1]
        result.DMA[3][:,6] = numpy.array(self.Oct)[:,1,1,2]
        result.DMA[3][:,7] = numpy.array(self.Oct)[:,0,2,2]
        result.DMA[3][:,8] = numpy.array(self.Oct)[:,1,2,2]
        result.DMA[3][:,9] = numpy.array(self.Oct)[:,0,1,2]
        #
        if self.has_hexadecapoles:
           result.DMA[4][:, 0] = numpy.array(self.Hex)[:,0,0,0,0]
           result.DMA[4][:, 1] = numpy.array(self.Hex)[:,1,1,1,1]
           result.DMA[4][:, 2] = numpy.array(self.Hex)[:,2,2,2,2]
           result.DMA[4][:, 3] = numpy.array(self.Hex)[:,0,0,0,1]
           result.DMA[4][:, 4] = numpy.array(self.Hex)[:,0,0,0,2]
           result.DMA[4][:, 5] = numpy.array(self.Hex)[:,1,1,1,0]
           result.DMA[4][:, 6] = numpy.array(self.Hex)[:,1,1,1,2]
           result.DMA[4][:, 7] = numpy.array(self.Hex)[:,2,2,2,0]
           result.DMA[4][:, 8] = numpy.array(self.Hex)[:,2,2,2,1]
           result.DMA[4][:, 9] = numpy.array(self.Hex)[:,0,0,1,1]
           result.DMA[4][:,10] = numpy.array(self.Hex)[:,0,0,2,2]
           result.DMA[4][:,11] = numpy.array(self.Hex)[:,1,1,2,2]
           result.DMA[4][:,12] = numpy.array(self.Hex)[:,0,0,1,2]
           result.DMA[4][:,13] = numpy.array(self.Hex)[:,1,1,0,2]
           result.DMA[4][:,14] = numpy.array(self.Hex)[:,2,2,0,1]

        # finally change the origins from numpy.zeros to origins (for C(A+B)MMs and LMTP)
        # tip: in the future CAMM should be rewritten in a vectorized manner as in CABMMs and LMTPs,
        #      which are computed at the origin of global coordinate system [0,0,0] and then recentered 
        #      to atomic or LMO centroid centers.
        #      Therefore, changing origins would be done also for CAMMs.
        if change_origins:
           result.MAKE_FULL()
           result.ChangeOrigin(new_origin_set=self.origin)
       
           # save the memorials after recenterization  
           result.MAKE_FULL()
           self.Dip = result.DMA_FULL[2]
           self.Quad= result.DMA_FULL[3]
           self.Oct = result.DMA_FULL[4]
           if self.has_hexadecapoles:
              self.Hex = result.DMA_FULL[5] 

        result.set_structure(pos=self.pos)
        self.__dma_bin.append(result)
        
        return
        
    def MU(self):
        """calculate molecular dipole moment""" 
        
        Mu = numpy.array([0,0,0], dtype=numpy.float64)
        # nuclear contribution
        if not self.transition:
           for atom in self.molecule.atoms:
               Mu+= atom.atno*numpy.array(atom.pos())
        # electronic contribution
        #Mu[:] -= self.P * self.D[:]
        for a in xrange(self.K):
            for b in xrange(self.K):
                Mu[0]-= self.P[a,b] * self.D[0,a,b]
                Mu[1]-= self.P[a,b] * self.D[1,a,b]
                Mu[2]-= self.P[a,b] * self.D[2,a,b]
        return Mu
    
    def QUAD(self):
        """calculate molecular quadrupole moment"""
        
        Quad = numpy.zeros( (3,3), dtype=numpy.float64)
        # nuclear contribution
        if not self.transition:
           for atom in self.molecule.atoms:
               Z = atom.atno
               R = numpy.array(atom.pos())
               for i in [0,1,2]:
                   for j in [0,1,2]:
                       Quad[i,j] += Z*R[i]*R[j]
        # electronic contribution
        for a in xrange(self.K):
            for b in xrange(self.K):  
                for i in [0,1,2]:
                    for j in [0,1,2]:
                        Quad[i,j] -= self.P[a,b] * self.Q[i,j,a,b]
        return Quad

    def OCT(self):
        """calculate molecular octupole moment"""
        
        Oct = numpy.zeros( (3,3,3), dtype=numpy.float64)
        # nuclear contribution
        if not self.transition:
           for atom in self.molecule.atoms:
               Z = atom.atno
               R = numpy.array(atom.pos())
               for i in [0,1,2]:
                   for j in [0,1,2]:
                       for k in [0,1,2]:
                           Oct[i,j,k] += Z*R[i]*R[j]*R[k]
        # electronic contribution
        for a in xrange(self.K):
            for b in xrange(self.K):  
                for i in [0,1,2]:
                    for j in [0,1,2]:
                        for k in [0,1,2]:
                            Oct[i,j,k] -= self.P[a,b] * self.O[i,j,k,a,b]
        #self.makeTracelessOCT(Oct)
        return Oct 

    def HEX(self):
        """calculate molecular hexadecapole moment"""
        
        Hex = numpy.zeros( (3,3,3,3), dtype=numpy.float64)
        # nuclear contribution
        if not self.transition:
           for atom in self.molecule.atoms:
               Z = atom.atno
               R = numpy.array(atom.pos())
               for i in [0,1,2]:
                   for j in [0,1,2]:
                       for k in [0,1,2]:
                           for l in [0,1,2]:
                               Hex[i,j,k,l] += Z*R[i]*R[j]*R[k]*R[l]
        # electronic contribution
        for a in xrange(self.K):
            for b in xrange(self.K):  
                for i in [0,1,2]:
                    for j in [0,1,2]:
                        for k in [0,1,2]:
                            for l in [0,1,2]:
                                Hex[i,j,k,l] -= self.P[a,b] * self.H[i,j,k,l,a,b]
        #self.makeTracelessHEX(Hex)
        return Hex


    def makeTracelessOCT(self,O):
        W = O.copy()
        O *= (5./2.)
        for i in [0,1,2]:
            for j in [0,1,2]:
                for k in [0,1,2]:
                    Wt = trace(W)
                    if   i==j and i==k:
                         O[i,j,k]-= (1./2.) * (Wt[i] + Wt[j] + Wt[k])
                    elif i==j and i!=k:
                         O[i,j,k]-= (1./2.) * Wt[k]
                    elif j!=k and i==k:
                         O[i,j,k]-= (1./2.) * Wt[j]
                    elif i!=j and j==k:
                         O[i,j,k]-= (1./2.) * Wt[i]
        return  


    def makeTracelessHEX(self,H):
        """
Formula: 

H_{ijkl}^{Traceless} = \frac{35}{8} H_{ijkl} 
                     - \frac{5}{8} ( \delta_{ij} H_{00kl} + \delta_{ik} H_{00jl} + \delta_{il} H_{00jk} + 
                                     \delta_{jk} H_{00il} + \delta_{jl} H_{00ik} + \delta_{kl} H_{00ij} )
                     + \frac{1}{8} Tr^2(H) (\delta_{ij} \delta_{kl} + \delta_{ik} \delta_{jl} + \delta_{il} \delta_{jk})

where 
 
  H_{ijkl} refers to primitive hexadecapole moment tensor element


       def.
  Tr^2 ==== XXXX + YYYY + ZZZZ + 2 (XXYY + XXZZ + YYZZ)


"""
        W = H.copy()
        H *= (35./8.)
        Wt = W.trace()
        Wtt= W.trace().trace()
        c = (5./8.)
        p = (1./8.)
        d = numpy.identity(3,numpy.float64)
        for i in [0,1,2]:
            for j in [0,1,2]:
                for k in [0,1,2]:
                    for l in [0,1,2]:
                        H[i,j,k,l] -= c * (d[i,j] * Wt[k,l] + d[i,k] * Wt[j,l] + d[i,l] * Wt[j,k] +
                                           d[j,k] * Wt[i,l] + d[j,l] * Wt[i,k] + d[k,l] * Wt[i,j] ) \
                                   -  p * Wtt * (d[i,j] * d[k,l] + d[i,k] * d[j,l] + d[i,l] * d[j,k] )
#        for i in [0,1,2]:
#            for j in [0,1,2]:
#                for k in [0,1,2]:
#                    for l in [0,1,2]:
#                        if i==j==k==l:                 # AAAA
#                           H -= (15./4.) * Wt[i,i]       - 3.0 * Wtt
#                        elif i==j and k==l and i!=k:   # AABB
#                           H -= c * (Wt[i,i] + Wt[k,k])  - Wtt
#                        elif i==k and j==l and i!=j:   # ABAB
#                           H -= c * (Wt[i,i] + Wt[j,j])  - Wtt
#                        elif i==l and j==k and i!=j:   # ABBA
#                           H -= c * (Wt[i,i] + Wt[j,j])  - Wtt
#                        elif i==j and j!=k and k!=l:   # AABC
#                           H -= c * Wt[k,l]
#                        elif i!=j and j!=k and k==l:   # ABCC
#                           H -= c * Wt[i,j]
#                        elif i!=j and j==k and k!=l:   # ABBC
#                           H -= c * Wt[i,l]
#                        elif i!=j and i==k and k!=l:   # ABAC
#                           H -= c * Wt[j,l]
#                        elif i!=j and j!=k and j==l:   # ABCB
#                           H -= c * Wt[i,k]
#                        elif i==l and i!=j and j!=k:   # ABCA
#                           H -= c * Wt[j,k]
        return  

    
    def makeTracelessCAMMs(self):
        """turns ordinary cartesian CAMMs into traceless cartesian CAMMs"""
        
        # Quadrupole moment
        for Q in self.Quad:
            t=Q.trace()
            Q *= (3./2.)
            for i in [0,1,2]: Q[i,i]-=(1./2.)*t
        # Octapole moment
        for O in self.Oct:
            W = O.copy()
            Wt= W.trace()
            O *= (5./2.)
            for i in [0,1,2]:
                for j in [0,1,2]:
                    for k in [0,1,2]:
                        if   i==j and i==k:
                             O[i,j,k]-= (1./2.) * (Wt[i] + Wt[j] + Wt[k])
                        elif i==j and i!=k:
                             O[i,j,k]-= (1./2.) * Wt[k]
                        elif j!=k and i==k:
                             O[i,j,k]-= (1./2.) * Wt[j]
                        elif i!=j and j==k:
                             O[i,j,k]-= (1./2.) * Wt[i]
        # Hexadecapole moment
        if self.has_hexadecapoles:
           #raise NotImplementedError, 'traceless hexadecapoles are not implemented yet!'
           for H in self.Hex:
               W = H.copy()
               H *= (35./8.)
               Wt = W.trace()
               Wtt= W.trace().trace()
               c = (5./8.)
               p = (1./8.)                                                                                    
               d = numpy.identity(3,numpy.float64)
               for i in [0,1,2]:
                   for j in [0,1,2]:
                       for k in [0,1,2]:
                           for l in [0,1,2]:
                               H[i,j,k,l] -= c * (d[i,j] * Wt[k,l] + d[i,k] * Wt[j,l] + d[i,l] * Wt[j,k] +
                                                  d[j,k] * Wt[i,l] + d[j,l] * Wt[i,k] + d[k,l] * Wt[i,j] ) \
                                          -  p * Wtt * (d[i,j] * d[k,l] + d[i,k] * d[j,l] + d[i,l] * d[j,k] )

        self.clock.actualize('transformation to traceless tensors')
             
    
    def ReturnCAMMs(self):
        """returns CAMMs"""
       
        if not self.has_hexadecapoles: 
           return self.RArray, self.Mon, self.Dip, self.Quad, self.Oct 
        else:
           return self.RArray, self.Mon, self.Dip, self.Quad, self.Oct, self.Hex


    def __printMMMs__(self):
        """
  .................................................................
  :                 MOLECULAR MULTIPOLE MOMENTS                   :
  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
"""

        muu   = self.DipoleMoment
        mud   = muu/0.39343029
        mu_av = numpy.sqrt(muu[0]**2 + muu[1]**2 + muu[2]**2)
        mud_av= mu_av/0.39343029

        _mu_ = "%12s\n"                           % ('[A.U.]'.rjust(12))
        _mu_+= "%12s %12s %12s %12s\n"            % ('X'.rjust(12),'Y'.rjust(12),'Z'.rjust(12),'ABS'.rjust(12)) 
        _mu_+= "%12.6f %12.6f %12.6f %12.6f\n\n"  % (muu[0], muu[1], muu[2], mu_av )
        _mu_+= "%12s\n"                           % ('[DEBYE]'.rjust(12))
        _mu_+= "%12s %12s %12s %12s\n"            % ('X'.rjust(12),'Y'.rjust(12),'Z'.rjust(12),'ABS'.rjust(12)) 
        _mu_+= "%12.6f %12.6f %12.6f %12.6f\n\n"  % (mud[0], mud[1], mud[2], mud_av )
 
        quadd  = self.QuarupoleMoment

        _quad_ = "%12s %12s %12s\n"               % ('XX'.rjust(12),'XY'.rjust(12),'XZ'.rjust(12)) 
        _quad_+= "%12.6f %12.6f %12.6f\n"         % (quadd[0][0], quadd[0][1], quadd[0][2]) 
        _quad_+= "%12s %12s %12s\n"               % ('YX'.rjust(12),'YY'.rjust(12),'YZ'.rjust(12)) 
        _quad_+= "%12.6f %12.6f %12.6f\n"         % (quadd[1][0], quadd[1][1], quadd[1][2]) 
        _quad_+= "%12s %12s %12s\n"               % ('ZX'.rjust(12),'ZY'.rjust(12),'ZZ'.rjust(12))
        _quad_+= "%12.6f %12.6f %12.6f\n\n"       % (quadd[2][0], quadd[2][1], quadd[2][2])

        octt = self.OctupoleMoment
     
        _oct_ = "%12s %12s %12s      \n"          % ('XXX'.rjust(12),'YYY'.rjust(12),'ZZZ'.rjust(12))   
        _oct_+= "%12.6f %12.6f %12.6f\n"          % (octt[0][0][0], octt[1][1][1], octt[2][2][2]) 
        _oct_+= "%12s %12s %12s      \n"          % ('XXY'.rjust(12),'XXZ'.rjust(12),'YYZ'.rjust(12))  
        _oct_+= "%12.6f %12.6f %12.6f\n"          % (octt[0][0][1], octt[0][0][2], octt[1][1][2]) 
        _oct_+= "%12s %12s %12s      \n"          % ('XZZ'.rjust(12),'XYY'.rjust(12),'YZZ'.rjust(12))  
        _oct_+= "%12.6f %12.6f %12.6f\n"          % (octt[0][2][2], octt[0][1][1], octt[1][2][2])
        _oct_+= "%12s                \n"          % ('XYZ'.rjust(12))  
        _oct_+= "%12.6f              \n\n"        % (octt[0][1][2])

        log = '\n'
        log+= " =========================== \n"
        log+= " MOLECULAR MULTIPOLE MOMENTS \n"
        log+= " =========================== \n\n"

        log+= "           ------------- \n"
        log+= "           DIPOLE MOMENT \n"
        log+= "           ------------- \n\n"
        log+= _mu_ 

        log+= "           ----------------- \n"
        log+= "           QUADRUPOLE MOMENT \n"
        log+= "           ----------------- \n\n"
        log+= _quad_

        log+= "           --------------- \n"
        log+= "           OCTUPOLE MOMENT \n"
        log+= "           --------------- \n\n"
        log+= _oct_
 
        print log

    def __printCAMMs__(self):   # old printout form
        """
  .................................................................
  :                 W. A. SOKALSKI & R. A. POIRIER                :
  :                 ******************************                :
  :        " CUMULATIVE ATOMIC MULTIPOLE REPRESENTATION           :
  :             OF THE MOLECULAR CHARGE DISTRIBUTION              :
  :                 AND ITS BASIS SET DEPENDENCE "                :
  :         CHEMICAL PHYSICS LETTERS, VOL. 98, NO. 1, 1983        :
  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
"""
        log = ""
        log+= '\n'
        log+= " CUMULATIVE CHARGES\n"
        log+= " ------------------\n"
        log+= " %s %s\n"                         % ('ATOM'.rjust(10), 'q'.rjust(12))
        for i in range(self.N_at):
            log+= " %10d %12.8f\n"               % ( i, self.Mon[i] )  
        log+= " \n"   
 
        log+= " CUMULATIVE DIPOLES\n"
        log+= " ------------------\n"
        log+= " %s %s %s %s\n"                   % ('ATOM'.rjust(10), 'x'.rjust(12), 'y'.rjust(12), 'z'.rjust(12))
        for i in range(self.N_at):
            log+= " %10d %12.8f %12.8f %12.8f\n" % ( i, self.Dip[i][0], self.Dip[i][1], self.Dip[i][2] )
        log+= " \n" 
 
        log+= " CUMULATIVE QUADRUPOLES\n"
        log+= " ----------------------\n"
        log+= " %s %s %s %s %s %s %s %s %s %s\n" % ('ATOM'.rjust(10), 'xx'.rjust(12), 'xy'.rjust(12), 'xz'.rjust(12),
                                                                      'yx'.rjust(12), 'yy'.rjust(12), 'yz'.rjust(12),
                                                                      'zx'.rjust(12), 'zy'.rjust(12), 'zz'.rjust(12))
        for i in range(self.N_at):
            log+= " %10d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n" % ( i, self.Quad[i][0][0], self.Quad[i][0][1], self.Quad[i][0][2], 
                                                                                                  self.Quad[i][1][0], self.Quad[i][1][1], self.Quad[i][1][2],
                                                                                                  self.Quad[i][2][0], self.Quad[i][2][1], self.Quad[i][2][2])
        log+= " \n" 
 
        log+= " CUMULATIVE OCTUPOLES\n"
        log+= " --------------------\n"
        log+= " %s %s %s %s %s %s %s %s %s %s\n" % ('ATOM'.rjust(10), 'xxx'.rjust(12), 'xxy'.rjust(12), 'xxz'.rjust(12),
                                                                      'xyx'.rjust(12), 'xyy'.rjust(12), 'xyz'.rjust(12),
                                                                      'xzx'.rjust(12), 'xzy'.rjust(12), 'xzz'.rjust(12))
        for i in range(self.N_at):
            log+= " %10d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n" % ( i, self.Oct[i][0][0][0], self.Oct[i][0][0][1], self.Oct[i][0][0][2], 
                                                                                                  self.Oct[i][0][1][0], self.Oct[i][0][1][1], self.Oct[i][0][1][2],
                                                                                                  self.Oct[i][0][2][0], self.Oct[i][0][2][1], self.Oct[i][0][2][2])
        log+= " \n" 

        log+= " %s %s %s %s %s %s %s %s %s %s\n" % ('ATOM'.rjust(10), 'yxx'.rjust(12), 'yxy'.rjust(12), 'yxz'.rjust(12),
                                                                      'yyx'.rjust(12), 'yyy'.rjust(12), 'yyz'.rjust(12),
                                                                      'yzx'.rjust(12), 'yzy'.rjust(12), 'yzz'.rjust(12))
        for i in range(self.N_at):
            log+= " %10d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n" % ( i, self.Oct[i][1][0][0], self.Oct[i][1][0][1], self.Oct[i][1][0][2], 
                                                                                                  self.Oct[i][1][1][0], self.Oct[i][1][1][1], self.Oct[i][1][1][2],
                                                                                                  self.Oct[i][1][2][0], self.Oct[i][1][2][1], self.Oct[i][1][2][2])
        log+= " \n" 

        log+= " %s %s %s %s %s %s %s %s %s %s\n" % ('ATOM'.rjust(10), 'zxx'.rjust(12), 'zxy'.rjust(12), 'zxz'.rjust(12),
                                                                      'zyx'.rjust(12), 'zyy'.rjust(12), 'zyz'.rjust(12),
                                                                      'zzx'.rjust(12), 'zzy'.rjust(12), 'zzz'.rjust(12))
        for i in range(self.N_at):
            log+= " %10d %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f %12.8f\n" % ( i, self.Oct[i][2][0][0], self.Oct[i][2][0][1], self.Oct[i][2][0][2], 
                                                                                                  self.Oct[i][2][1][0], self.Oct[i][2][1][1], self.Oct[i][2][1][2],
                                                                                                  self.Oct[i][2][2][0], self.Oct[i][2][2][1], self.Oct[i][2][2][2])
        log+= " \n" 

  
        print log

    def __repr__(self):
        """new printout form"""
        log = ''
        for dma in self.__dma_bin:
            log+=str(dma)
        return log
