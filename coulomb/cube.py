#*-* coding: utf-8 *-*
# --------------------------------- #
#               CUBE                #
# --------------------------------- #

__all__=['CUBFLE','QMMap']

from numpy import *
from libbbg.utilities import Read_xyz_file, Vr_dma, Energy_density, SVDSuperimposer
from libbbg.units import *
from PyQuante.GridPoint import GridPoint
import sys

class CUBFLE:
      """
 ---------------------------------------------------------------------------------------------
 Represents Gaussian CUBE file

                                                                Last Revision: 14 Dec 2014
 ---------------------------------------------------------------------------------------------

 Usage:

  a = CUBFLE(file_name)

  dim   = a.get_dimensions()
  coord = a.get_coord()
  data  = a.get_data()
  sp    = a.get_spacings()
  orig  = a.get_origin()

  a.rotate(rot)
  a.translate(transl)
  a.sup(xyz, supl=None)
  a.write(name='happy.cube', misc='comment')

  int = a.integrate(rule='simple')

  print a

 Notes:
  
  all quantities are stored and returned in A.U.! 

 ---------------------------------------------------------------------------------------------
                               Author: Bartosz Błasiak, email: <blasiak.bartosz@gmail.com>
 ---------------------------------------------------------------------------------------------
"""
      
      def __init__(self, file, discalloc=False):
          self.file = file
          (self.dimensions, self.spacings, self.origin , 
           self.N_atoms,   self.coord)  = self.read_basic(file)
          if not discalloc: self.data   = self.read_cube(file)
          self.n_x, self.n_y, self.n_z  = self.dimensions

      def get_dimensions(self):
          return self.dimensions
      def get_coord(self):
          return self.coord
      def get_data(self):
          return self.data
      def get_spacings(self): 
          return self.spacings
      def get_origin(self):
          return self.origin

      def integrate(self, rule='simple'):
          """integrate the cube file content"""
          s = linalg.eig(self.spacings)[0]
          Int = s[0]*s[1]*s[2]
          Io = 0.0
          #
          if rule=='simple': Io = self.data.sum()
          else: raise Exception('No rule %s implemented' % rule)
          return Io * Int 

      def sup(self, xyz, suplist=None, rotran=None):
          """superomposes the cube file content"""
          if rotran is None:
             if len(xyz)!=1:
                s = SVDSuperimposer()
                if suplist is None: s.set(xyz,self.__pos)
                else:               s.set(xyz[suplist],self.__pos[suplist])
                s.run()
                rms         = s.get_rms()
                rot, transl = s.get_rotran()
             else:  # 1-atomic case
                rot = identity(3)
                transl = xyz - self.__pos 
                rms = 0.0
          else:
             rot, transl = rotran
             rms = 0.0
          self.rotate(rot)
          self.translate(transl)
          return rms

      def rotate(self, rot):
          """rotates the cube file content"""
          self.origin = dot(self.origin, rot)
          self.spacings = dot(self.spacings, rot)
          self.coord[:,2:] = dot(self.coord[:,2:], rot)
          return

      def translate(self, transl):
          """translates the cube file content"""
          self.origin+= transl
          self.coord[:,2:] += transl
          return

      def __repr__(self):
          """Pring me!"""
          hline = ' '+31*'-'+'\n'
          log = '\n'
          log+= hline
          log+= ' g09 Cube file\n'
          log+= hline 
          # file name
          log+= ' %10s %20s\n' % ('Name   :'.ljust(10), self.file.rjust(20))
          # atoms
          log+= ' %10s %20d\n' % ('N_atoms:'.ljust(10), self.N_atoms)
          # data size
          log+= ' %10s %20d\n' % ('Size   :'.ljust(10), len(self.data))
          log+= hline
          return str(log)

      def read_inst(self):
          """?"""
          pass
        
      def read_basic(self,file):
          """reads .cube file and parses info"""
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

          # cube dimensions and shape
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
 
          dimensions = array([n_x,n_y,n_z])
          spacings  = array([vox_x,vox_y,vox_z])

          return dimensions, spacings, origin, N_atoms, coord

      def read_cube(self,file):
          """reads .cube files and parses info"""
          File = open(file)
          for i in range(6+self.N_atoms): File.readline()
          # function values
          data = []
          line = File.readline() 
          while line:
                data+= line.split()
                line = File.readline()
          data = array(data,dtype=float64)
          return data

      def write(self, name, misc='No comments on this file'):
          """writes the cube file"""
          self._write_cube(name, misc, self.N_atoms, self.origin, 
                           self.dimensions, self.spacings, self.coord, self.data)
          return

      def write_ext(self, name, N_atoms, origin, dimensions, spacings, coord, data, misc='No comments on this file'):
          """writes the cube file by reading external data"""
          self._write_cube(name, misc, N_atoms, origin, 
                           dimensions, spacings, coord, data)
          return

      # PROTECTED

      def _write_cube(self, name, misc, N_atoms, origin, dimensions, spacings, coord, data):
          """Write cube file from external source"""
          out = open(name,'w')
          print >> out, " File: %s\n %s" % (name, misc)
          print >> out, "%5d" % N_atoms, 
          print >> out, 3*"%11.6f " % tuple(origin)
          for i in range(3):
              print >> out, "%5d" % dimensions[i],
              print >> out, 3*"%11.6f " % tuple(spacings[i])
          for i in range(N_atoms):
              print >> out, "%5d" % coord[i][0],
              for c in range(len(coord[0])-1):
                  print >> out, "%11.6f" % coord[i][c+1],
              print >> out

          k = 1
          d = dimensions[2]
          log = ""
          for i in range(len(data)):
              if   (k%6!=0 and k%d!=0):
                    k+=1
                    log+="%13.5E" % data[i]
              elif (k%6==0 and k%d!=0):
                    log+="%13.5E\n" % data[i]
                    k+=1
              elif (k%6!=0 and k%d==0):
                    log+="%13.5E\n" % data[i]
                    k=1
              elif (k%6==0 and k%d==0):
                    log+="%13.5E\n" % data[i]
                    k=1
          print >> out, log
          out.close()
          return
        
class QMMap(CUBFLE):
    """
 ---------------------------------------------------------------------------------------------
 Molecular Quantum Map

                                                                 Last Revision: 14 Dec 2014
 ---------------------------------------------------------------------------------------------

 I. DESCRIPTION.

    Represents a grid map of some molecular space distribution.                              
    The data can be generated:
                                                                                             
      o directly from the density matrix and basis set, or 
      o from DMA distribution. 
                                                                                             
    The result is written to the g09 cube file.
                                                                                             
    The available types of distributions are shown in Table I.
    
    Table I. Types of distribution available in QMMap
    --------------------------------------------------------------------
      Type of cube                                    WFN   DMA
    --------------------------------------------------------------------
      distribution density                            YES   NO 
      electrostatic (generalized) potential           YES   YES
      electrostatic (generalized) energy density      NO    YES
    --------------------------------------------------------------------
      * YES - implemented ; NO - not implemented yet
    --------------------------------------------------------------------
    
    If the distributions are generated from DMA, points within atomic vdW radii are 
    discarded (zero values in these points are written because no damping implemented yet). 

 II. USAGE.

     a = QMMap( coord, dmat, bfs, dma, dimensions=(60,60,60), padding=(10.0,10.0,10.0))
     a.eval('dens')
     a.write('map.cube', misc='My first QMMap')

     print a.npoints
     print a.h  # increment in the grid

 III. ARGUMENTS.
     
     coord      - list of atomic symbols and coordinates (in Bohr) of the form:     
                  [ list_1, list_2, ..., list_n ]
                  where list_n = ['Na', 0.2, 0.3, 0.4] for example
                  One can generate after opening xyz file by libbbg.utilities.QMFile:
                  coord = QMFile_instance.get_coord_list()
     dmat       - ndarray of size nbs x nbs
     bfs        - PyQuante.BasisSet object of length nbs
                  It can be easily generated e.g. by QMFile:
                  bfs = QMFile_instance.get_bfs()
     dma        - libbbg.dma.DMA object
     dimensions - dimensions of the resulting cube (default is 60 x 60 x 60 points)
     padding    - padding (in Bohr) of the cube in x, y and z directions. The size 
                  of the cube will be extended by +/- u_max/u_min (u=x,y,z) where
                  u_min/max are the minimal/maximal vaules of atomic coordinates.
     npts       - uniform grid density of npts points per bohr. Use 3 for coarse, 5 for 
                  medium and 10 for dense grids.

 IV. NOTES.

     o This routine produces cubic grids. The increment in each direction is
       assumed to be constant and is calculated as:
       h = ( xyz_max - xyz_min) / array(dimensions)
       where xyz_min/max are the terminal coordinates after adding padding vectors
       Therefore the 'standard orientation' of molecule is recommended for the compact
       cubes (i.e., the main axes of symmetry being colinear with global coordinate axes).

     o Instead of defining the grid by providing dimensions vector, one may
       define the uniform grid density on npts point per bohr. In such a case the dimensions
       are calculated accordingly.

  V. SYNOPSIS.

     None.

 VI. SEE ALSO: 

     libbbg.utilities.QMFile
     PyQuante.Molecule
     libbbg.dma.DMA

 ---------------------------------------------------------------------------------------------
                               Author: Bartosz Błasiak, email: <blasiak.bartosz@gmail.com>
 ---------------------------------------------------------------------------------------------
"""
    
    def __init__(self, coord_list, dmat=None, bfs=None, dma=None, 
                       dimensions=(60,60,60), padding=(4.0,4.0,4.0),npts=None):
        self.file = 'None'
        # check the type of calculation
        if (dmat is None and dma is None): 
            raise TypeError(" You must specify either WFN or DMA to proceed! Termination.")
        if dmat is not None:
            self.calc_type = 'wfn'
        else:
            self.calc_type = 'dma'
        # gather the memorials 
        self.coord = coord_list
        self.N_atoms = len(coord_list)
        # replace atomic symbols by atomicnumbers and insert the zero column
        for i in range(self.N_atoms):
            self.coord[i][0] = Atom(self.coord[i][0]).atno
            self.coord[i].insert(1, 0.00000)
        self.coord = array(self.coord)
        self.dimensions=dimensions
        self.npoints = prod(self.dimensions)
        self.padding=array(padding,dtype=float64)
        self.npts = npts
        self.dmat = dmat
        self.bfs = bfs
        # prepare the dma object if any
        if dma is not None:
           self.dma = dma.copy()
           self.dma.MAKE_FULL()     
           self.dma.MakeTraceless()
        else: self.dma = None
        #
        self.__CalculateOrigin() ### ---> origin, xmin,xmax etc...
        
               
    def eval(self, which='potential'):
        """
Evaluates the space-distributed property
Usage:
 eval(which=<type>)
where <type> = 'potential', 'density' or 'energydensity'
Default: 'potential'

This method accumulates CUBFLE.data memorial
"""
        # DMA-based calculation
        if self.calc_type == 'dma':
           if which.lower()[:3] == 'pot':
              self.property = 'electrostatic potential'                              
              for i in range(self.npoints):
                  r = self.__point(i)
                  allowed = True
                  for atom in self.coord:
                      if self.__vdW_sphere(Atom(atom[0]).radius,r,array(atom[1:])): 
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
                  #    if self.__vdW_sphere(Atom(atom[0]).radius,r,array(atom[1:])):
                   #      allowed = False
                  #       continue
                  if allowed:
                     self.data[i] = Energy_density(self.dma,r,is_full=True)

        # WFN-based calculation
        elif self.calc_type == 'wfn':
           if which.lower()[:3] == 'den':
              self.property = 'distribution density'
              for i in range(self.npoints):
                  x,y,z = self.__point(i)
                  point = GridPoint(x,y,z)
                  point.set_bf_amps(self.bfs)
                  point.setdens(self.dmat)
                  self.data[i] = point.dens()

           elif which.lower()[:3] == 'pot':
             self.property = 'electrostatic potential'
             raise NotImplementedError
        
    # PRIVATE

    def __vdW_sphere(self,vdW_radius,R,Ri):
        """Checks if point lies within atom's vdW sphere. 
If '0' - the point lies outside atom's vdW sphere."""
        r2 = sum((R-Ri)**2,axis=0)
        if r2>vdW_radius**2: return 0
        else:                return 2

    def __CalculateOrigin(self):
        """Setup the grid data"""
        n=len(self.coord)
        d=zeros((n,3),dtype=float64)
        for i in range(n):
            d[i] = asarray(self.coord[i][2:])
        self.xyz_min = ndarray.min(d,axis=0)
        self.xyz_max = ndarray.max(d,axis=0)
        
        self.xyz_min -= self.padding
        self.xyz_max += self.padding
        
        self.origin = self.xyz_min
        
        #change dimensions for a regular grid of npts points per bohr
        self.box_size = self.xyz_max - self.xyz_min
        if self.npts is not None:
            self.h = ones(3)/self.npts
            self.dimensions = rint(self.box_size/self.h)
            self.dimensions = self.dimensions.astype(int)
            self.npoints = prod(self.dimensions)
        
        self.h = self.box_size/array(self.dimensions)
        self.spacings = diag(self.h)
        print self.h
        print self.spacings
        print self.box_size
        print self.dimensions
        
        self.data = zeros(self.dimensions[0]*self.dimensions[1]*self.dimensions[2],float64)
        
    def __point(self,i):
        """Generate the coordinates of (i+1)-th point in the Density Cube"""
        nx,ny,nz=self.dimensions
        ix = i/(ny*nz)
        iy = (i-ix*ny*nz)/nz
        iz = i-ix*ny*nz-iy*nz
        r = self.xyz_min + ix*self.spacings[0] + iy*self.spacings[1] + iz*self.spacings[2]
        return r

