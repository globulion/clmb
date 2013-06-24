# ------------------------------------------- #
#         MULTIPOLE DISTRIBUTION MODULE       #
# ------------------------------------------- #

from run  import *
from utilities2 import array_outer_product,    \
                       array_outer_product_1_2,\
                       array_outer_product_2_1

__all__=['MULTIP']

class MULTIP(RUN):
    """contains useful procedures for computing:
       - Molecular Multipole Moments (MMMs)
       - Cumulative Atomic Multipole Moments (CAMMs)"""
       
    def __init__(self, molecule, basis, method,matrix=None,multInts=None,transition=False):
        RUN.__init__(self, molecule, basis, method,matrix,multInts)
        # LIST1 - the list of atoms in the order of basis functions used, 
        # e.g.: for h2o molecule with atoms: 8,1,1 and STO-3G basis 
        # LIST1 = 0    0    0    0    0    1    2
        #         |    |    |    |    |    |    |
        #        1s   2s   2px  2py  2pz  1s   1s
        self.LIST1 = self.bfs.LIST1 
        self.transition = transition
        self.name=molecule.name

    def mmms(self):
        """calculate molecular multipole moments"""
        
        self.DipoleMoment = self.MU()
        self.QuarupoleMoment = self.QUAD()
        self.OctupoleMoment  = self.OCT()
        #self.HexadecapoleMoment = self.HEX()
        self.clock.actualize("calculation of MMMs")

    def camms(self):
        """calculate camms"""
        Mon      = []
        Dip      = []
        Quad     = []
        Oct      = []
        n = 0 # atom number
        for atom in self.molecule.atoms:
            if self.transition: q  = 0
            else:               q=atom.atno 
            M  = zeros( 3,dtype=float64)
            Q  = zeros((3,3),dtype=float64)
            O  = zeros((3,3,3),dtype=float64)
            R = array(atom.pos()) 
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
                       #R2D= array_outer_product_2_1(R2,self.D[:,I,J])
                       #RQ = array_outer_product_1_2(R,self.Q[:,:,I,J])
                       #O += self.P[I,J] * (\
                       #                   array_outer_product_1_2(R,outer(R,R)) * self.S[I,J] \
                       #                 - R2D - transpose(R2D,axes=(0,2,1)) - transpose(R2D,axes=(2,1,0))\
                       #                 + RQ  + transpose(RQ ,axes=(1,2,0)) + transpose(RQ ,axes=(2,0,1))\
                       #                 - self.O[:,:,:,I,J] 
                       #                   )  


            Mon.append(q)
            Dip.append(M)
            Quad.append(Q)
            Oct.append(O)
            n+=1     
        # save CAMMS 
        self.Mon   = Mon
        self.Dip   = Dip
        self.Quad  = Quad
        self.Oct   = Oct
 
        # clock measure
        self.clock.actualize('computing CAMMs')

    def MU(self):
        """calculate dipole moment""" 
        
        Mu = array([0,0,0],dtype=float64)
        # nuclear contribution
        if not self.transition:
           for atom in self.molecule.atoms:
               Mu+= atom.atno*array(atom.pos())
        # electronic contribution
        #Mu[:] -= self.P * self.D[:]
        for a in xrange(self.K):
            for b in xrange(self.K):
                Mu[0]-= self.P[a,b] * self.D[0,a,b]
                Mu[1]-= self.P[a,b] * self.D[1,a,b]
                Mu[2]-= self.P[a,b] * self.D[2,a,b]
        return Mu
    
    def QUAD(self):
        """calculate quadrupole moment"""
        
        Quad = zeros( (3,3) ,dtype=float64)
        # nuclear contribution
        if not self.transition:
           for atom in self.molecule.atoms:
               Z = atom.atno
               R = array(atom.pos())
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
        """calculate octupole moment"""
        
        Oct = zeros( (3,3,3) ,dtype=float64)
        # nuclear contribution
        if not self.transition:
           for atom in self.molecule.atoms:
               Z = atom.atno
               R = array(atom.pos())
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

    def makeTracelessOCT(self,O):
        W= zeros((3,3,3),dtype=float)
        W[:]= O
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
    
    def makeTracelessCAMMs(self):
        """turns ordinary cartesian CAMMs into traceless cartesian CAMMs"""
        
        # Quadrupole moment
        for Q in self.Quad:
            t=trace(Q)
            Q *= (3./2.)
            for i in [0,1,2]: Q[i,i]-=(1./2.)*t
        # Octapole moment
        for O in self.Oct:
            W= zeros((3,3,3),dtype=float)
            W[:]= O
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
        self.clock.actualize('transformation to traceless tensors')
             
    
    def ReturnCAMMs(self):
        """returns CAMMs"""
        
        return self.RArray, self.Mon, self.Dip, self.Quad, self.Oct 


    def __printMMMs__(self):
        """
  .................................................................
  :                 MOLECULAR MULTIPOLE MOMENTS                   :
  '''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
"""

        muu   = self.DipoleMoment
        mud   = muu/0.39343029
        mu_av = sqrt(muu[0]**2 + muu[1]**2 + muu[2]**2)
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

    def __printCAMMs__(self):
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
