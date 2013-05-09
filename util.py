#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  * --------------------------------------------------------- *
  *                   QUANTUM CHEMISTRY                       *
  *                  CALCULATION PACKAGE                      * 
  *                     C O U L O M B                         *
  * --------------------------------------------------------- *
  *  Author        :    Bartosz Błasiak & Robert W. Góra      *
  *  Affiliation   :    Wrocław University of Technology      *
  *  Modules used  :    PyQuante (modified), NumPy            *
  *  Libraries     :    LibInt (Valeev's 2-el integral code), *
  *                :    Lapack 1.4.0                          *
  * --------------------------------------------------------- *
  *  Options:                                                 *
  *    -h          :    print this help                       *
  *    -v          :    print description                     *
  *    [file]      :    load input file for calculations      *
  *                         - - -                             *
  * WARNING: YOU HAVE TO HAVE MODIFIED VERSION OF PYQUANTE    *
  *          INSTALLED. SEE THE MINI-MANUAL FOR DETAILS.      *
  *                                                           *
  *                                              @Globulion   *
  *                                                           *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
"""

# -------------------------------------------- #
#               U T I L I T I E S              #
# -------------------------------------------- #

from coulomb_head import *

__author__  = "Bartosz Błasiak & Robert W. Góra"
__version__ = "1.0"

__all__=['read_transition_dmatrix','TIMER']

def Error():
    """Error box""" 
    print
    print " ******************************************* "
    print " ***  ERROR:  Invalid option! :D         *** "
    print " ***     see the help box:               *** "
    print " ***          this_file.py -v            *** " 
    print " ******************************************* "
    print;exit()
 
 
# -------------------- #
#     COULOMB INFO     #
# -------------------- #

def Usage():
    """Print usage information and exit."""
    
    os.system('clear')
    print __doc__
    print " Machine epsilon is: ",finfo(float64).eps,"for float64 type\n"
    print;exit()

def Version():
    """Print version and other information"""
    
    print """
    COULOMB.py © 2012 version %s
    %s

    This simple script is designed to perform quantum chemi-   
    cal calculations of coulomb or pseudo-coulomb interaction  
    energy between molecules or molecular aggregates. The      
    methods implemented in this script are based on:           
    ---------------------------------------------------------  
     - DC          : Density Cube                              
     - ESP         : Charges from Electrostatic Potential      
     - CAMM        : Cumulative Atomic Multipole Moments       
     - EELEDS      : 1-order electrostatic energy term from 
                   : Hybrid Variational-Perturbational 
                   : Interaction Energy Decomposotion Scheme              
    ---------------------------------------------------------  
    It also can compute molecular multipole and CAMM distri-   
    butions as well as fit the charges according to ESP pro-
    cedure.

    The format of input file is described in mini-manual ava-  
    ilable in the package.                                     
   
    Usage:
    coulomb.py [file] - load input file (may contain molecule, units, methods,
                        density matrix and multipole ints)"
"""%(__version__,__author__)



def read_transition_dmatrix(matrix_querry='Alpha transition density to state',
                            filetype='gaussian',file=0,state=1):
    """reads transition density matrix from a .log file"""
    
    if filetype == 'gaussian':
       F = open(file)
       line = F.readline()
       while not 'basis functions,' in line: line = F.readline()
       nbf = int(line.split()[0])
       while not ((matrix_querry in line) and (str(state) in line)):
             line = F.readline()
       tdmat = zeros((nbf,nbf),dtype=float64)
       for x in range(nbf/5 + bool(nbf%5)): 
           line = F.readline()
           #indices = [ int(j) for j in line.split() ]
           n = int(line.split()[ 0]) - 1
           m = int(line.split()[-1]) #- 1
           for i in range(nbf):
               line = F.readline()
               numbers = [ float64(k.replace('D','E')) for k in line.split()[1:] ]
               tdmat[i,n:m] = numbers
       return tdmat
      
      
class TIMER:
    """process timing statistics utility class""" 
    
    def __init__(self):
        self.occurence = 'start'
        self.t0 = time()
        self.tp = self.t0
        self.occurence_list = []
        self.log = "\n"
        self.log+= " ------------------ \n"
        self.log+= " Timing information \n"
        self.log+= " ------------------ \n"
        self.log+= "\n" 
        self.total_time = 0
        
    def actualize(self,new_occurence):
        """actualizes new occurence in the clock history"""
        
        self.occurence = new_occurence
        self.measure()
        self.total_time = sum(self.occurence_list)
        
    def measure(self):
        """measures length of occurence"""
        
        self.tn = time()
        self.length = self.tn - self.tp
        self.occurence_list.append(self.length)
        self.log += " - %44s %30.5f sec\n" % (self.occurence.ljust(44),self.length)
        self.tp = self.tn
        
    def __print__(self):
        """prints itself"""
        
        suma = sum(self.occurence_list)
        self.log=self.log.split('\n')
        for i in range(len(self.log) - 6):
            self.log[i+5] += "(%4.1f%%)" % (self.occurence_list[i]/suma*100)
        for i in range(len(self.log)):
            print self.log[i]
        t = self.total_time
        print
        print    "  ================================================"
        print    "  TOTAL TIME:  %d days %d hours %d min %d sec " % ( t/86400,t/3600,t/60,int(t) )
        print;print
        
def PRINTL(M,list1="",list2=""):
    """ print helper 2 """
    d = 5
    L = len(transpose(M))

    if L % d == 0:
       n = L / d
    else:
       n = L / d + 1

    if list1 == '' and list2 == '':
       for b in range(n):
           try:
               m = M[:,(b*d):((b+1)*d)]
           except IndexError:
               m = M[:,(b*d):-1]

           for u in range(len(m)):
             for i in range(len(transpose(m))):
               v = "%12.5f" % m[u][i]
               print "%14s" % v.rjust(14),
             print
           print

    elif list2 == '':

       for b in range(n):
           try:
               l1 = list1[(b*d):((b+1)*d)]
               m  = M[:,(b*d):((b+1)*d)]
           except IndexError:
               l1 = list1[(b*d):-1]
               m = M[:,(b*d):-1]

           for i in range(len(l1)):
               t = "%s" % l1[i]
               print "%10s" % t.rjust(10),
           print
           for i in range(len(l1)):
               kk = '-'*8
               print "%s" % kk.rjust(10),
           print

           for u in range(len(m)):
             for i in range(len(transpose(m))):
               v = "%.6e" % m[u][i]
               print "%10s" % v.rjust(10),
             print
           print

    else:

       for b in range(n):
           try:
               l1 = list1[(b*d):((b+1)*d)]
               l2 = list2[(b*d):((b+1)*d)]
               m  =   M[:,(b*d):((b+1)*d)]
           except IndexError:
               l1 = list1[(b*d):-1]
               l2 = list2[(b*d):-1]
               m =    M[:,(b*d):-1]

           for i in range(len(l1)):
               t1 = "%s" % l1[i]
               print "%15s" % t1.rjust(15),
           print
           for i in range(len(l1)):
               t2 = "%s" % l2[i]
               print "%15s" % t2.rjust(15),
           print
           for i in range(len(l1)):
               kk = '-'*13
               print "%s" % kk.rjust(15),
           print

           for u in range(len(m)):
             for i in range(len(transpose(m))):
               v = "%.6e" % m[u][i]
               print "%15s" % v.rjust(15),
             print
           print
           
           
           
def Periodic(mendeleiev):
    '''Returns the mendeleiev table as a python list of tuples. Each cell
    contains either None or a tuple (symbol, atomic number), or a list of pairs
    for the cells * and **. Requires: "import re". Source: Gribouillis at
    www.daniweb.com - 2008 '''

    # L is a consecutive list of tuples ('Symbol', atomic number)
    L = [ (e,i+1) for (i,e) in enumerate( re.compile ("[A-Z][a-z]*").findall('''
    HHeLiBeBCNOFNeNaMgAlSiPSClArKCaScTiVCrMnFeCoNiCuZnGaGeAsSeBrKr
    RbSrYZrNbMoTcRuRhPdAgCdInSnSbTeIXeCsBaLaCePrNdPmSmEuGdTbDyHoEr
    TmYbLuHfTaWReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFm
    MdNoLrRfDbSgBhHsMtDsRgUubUutUuqUupUuhUusUuo'''))]

    # The following fills the void with nones and returns the list of lists
    mendeleiev = 0

    if mendeleiev:
        for i,j in ( (88,103), (56,71) ):
            L[i] = L[i:j]
            L[i+1:] = L[j:]
        for i,j in ( (12,10 ), (4,10), (1,16) ):
            L[i:i]=[None]*j 

        return [ L[18*i:18*(i+1)] for i in range (7) ]

    # Return a plain list of tuples
    else:
        return L

def Atomn(s,ptable):
    '''Returns the atomic number based on atomic symbol string
    ptable is a list of consecutive (symbol, atomic number) tuples.'''
    for n,a in enumerate(ptable):
        if a[0].lower().find(s.strip().lower()) !=-1 :
            return float(n+1)
