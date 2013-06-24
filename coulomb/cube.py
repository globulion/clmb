# --------------------------------- #
#               CUBE                #
# --------------------------------- #

__all__=['CUBFLE','PotentialMap']

from numpy import *
from utilities import Read_xyz_file, Vr_dma, Energy_density
from units import *

class CUBFLE:
      """represents Gaussian CUBE file"""
      
      def __init__(self,file,discalloc=False):
          self.file = file
          (self.dimension, self.spacings, self.origin , 
           self.N_atoms,   self.coord)  = self.read_basic(file)
          if not discalloc: self.data   = self.read_cube(file)
          self.n_x, self.n_y, self.n_z  = self.dimension
          
      def read_inst():
          """?"""
          pass
        
      def read_basic(self,file):
          """reads .cube files and parses info"""
          
          File = open(file)
          for i in range(2): File.readline()

          # first set-up lines
          L3      = File.readline()
          L4      = File.readline()
          L5      = File.readline()
          L6      = File.readline()

          # basic info
          N_atoms =   int(L3.split()[0])
          origin  = array(L3.split()[1:],dtype=float64)

          # cube dimension and shape
          n_x     =   int(L4.split()[0])
          vox_x   = array(L4.split()[1:],dtype=float64)
          n_y     =   int(L5.split()[0])
          vox_y   = array(L5.split()[1:],dtype=float64)
          n_z     =   int(L6.split()[0])
          vox_z   = array(L6.split()[1:],dtype=float64)

          # atomic coordinates 
          coord   = zeros( (N_atoms,5) )
          for i in range(N_atoms):
              line = array(File.readline().split(), dtype=float64)
              coord[i,:] = line

          return (n_x,n_y,n_z) , (vox_x,vox_y,vox_z) , origin, N_atoms, coord

      def read_cube(self,file):
          """reads .cube files and parses info"""
          
          File = open(file)
          for i in range(6): File.readline()

          # function values
          data = []
          line = File.readline() 
          while line:
                Line = line.split()
                for j in Line:
                    data.append(j)
                line = File.readline()

          data = array(data,dtype=float64)

          return data

      def write(self, name, N_atoms, origin, D, V, coord, Data):
          """writes the cube file"""
          
          out = open(name,'w')
          print >> out, " File blaaa\n blaaa!!!"
          print >> out, "%4d%12.6f%12.6f%12.6f" % (N_atoms, origin[0], origin[1], origin[2])
          for i in range(3):
              print >> out, "%4d%12.6f%12.6f%12.6f" % (D[i], V[i][0], V[i][1], V[i][2])
          for i in range(N_atoms):
              print >> out, "%4d" % Atom(coord[i][0]).atno,
              for c in range(len(coord[0])-1):
                  print >> out, "%11.6f" % coord[i][c+1],
              print >> out

          k = 1
          d = D[2]
          log = ""
          for i in range(len(Data)):
              if   (k%6!=0 and k%d!=0):
                    k+=1
                    log+="%12.5E " % Data[i]
              elif (k%6==0 and k%d!=0):
                    log+="%12.5E \n" % Data[i]
                    k+=1
              elif (k%6!=0 and k%d==0):
                    log+="%12.5E \n" % Data[i]
                    k=1
              elif (k%6==0 and k%d==0):
                    log+="%12.5E\n" % Data[i]
                    k=1
          print >> out, log
          out.close()

class PotentialMap(CUBFLE):
    """represents a potential map"""
    
    def __init__(self,xyzfile,dma,nxyz=(60,60,60),padding=(10,10,10)):
        self.coord = Read_xyz_file(xyzfile)
        self.npoints = nxyz[0]*nxyz[1]*nxyz[2]
        dma.MAKE_FULL()
        dma.MakeTraceless()
        self.dma = dma
        self.nxyz=nxyz
        self.padding=array(padding,dtype=float64)
        self.__CalculateOrigin() ### ---> origin, xmin,xmax etc...
        
    def __CalculateOrigin(self):
        n=len(self.coord)
        d=zeros((n,3),dtype=float64)
        for i in range(n):
            d[i] = asarray(self.coord[i][1:])
        self.xyz_min = ndarray.min(d,axis=0)
        self.xyz_max = ndarray.max(d,axis=0)
        
        self.xyz_min -= self.padding
        self.xyz_max += self.padding
        
        self.origin = self.xyz_min#zeros(3)#1./2. * (self.xyz_max + self.xyz_min)
        
        self.h = (self.xyz_max - self.xyz_min)/array(self.nxyz)
        self.vox = zeros((3,3),dtype=float64)
        self.vox = diag(self.h)
        
        self.data = zeros(self.nxyz[0]*self.nxyz[1]*self.nxyz[2])
        
    def __point(self,i):
        """"""
        nx,ny,nz=self.nxyz
        ix = i/(ny*nz)
        iy = (i-ix*ny*nz)/nz
        iz = i-ix*ny*nz-iy*nz
        r = self.xyz_min + ix*self.vox[0] + iy*self.vox[1] + iz*self.vox[2]
        return r
               
        
    def eval(self,pot=False):
        """evaluates eletrostatic potential for given molecule"""
        print "\n Evaluation of %s on %d points ..." % (str(where(pot,'electrostatic potential',
                                                                      'energy density')),
                                                                      self.npoints,)
        if pot:
           self.property = 'electrostatic potential'
           for i in range(self.npoints):
               r = self.__point(i)
               allowed = True
               for atom in self.coord:
                   if self.vdW_sphere(Atom(atom[0]).radius,r,array(atom[1:])): 
                      allowed = False
                      continue
               if allowed:
                  self.data[i] = Vr_dma(self.dma,r,is_full=True)
        else:
           self.property = 'electrostatic energy density'
           for i in range(self.npoints):
               r = self.__point(i)
               allowed = True
               #for atom in self.coord:
               #    if self.vdW_sphere(Atom(atom[0]).radius,r,array(atom[1:])):
                #      allowed = False
               #       continue
               if allowed:
                  self.data[i] = Energy_density(self.dma,r,is_full=True)

        
    def tofile(self,name):
        """"""
        self.write(name, len(self.coord), self.origin, self.nxyz, self.vox, self.coord, self.data)
        print "\n Cube file writthen to < %s > !" % name
        
    def vdW_sphere(self,vdW_radius,R,Ri):
        """checks if point lies within atom's vdW sphere."""
        """if '0' - the point lies outside atom's vdW sphere."""
       
        r2 = sum((R-Ri)**2,axis=0)
        if r2>vdW_radius**2: return 0
        else:                return 2
