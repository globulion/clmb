# --------------------------------------------- #
#               WORKING UTILITIES               #
# --------------------------------------------- #

from libbbg.units  import *
from libbbg.dma    import DMA 
from parser import PARSER
from multip import * 
from esp    import *
from eeleds import *
from eint   import * 
from numpy  import array

__all__=['DO']

class DO(PARSER):
    """contains main routine for COULOMB.py"""
    outname = None
    def __init__(self,file):
        # --- prepare commands to be invoked
        PARSER.__init__(self,file)
        # run
        #RUN(self.molecule,self.basis,self.method,matrix=self.dmat,multInts=self.multints)
        #run.run(self.molecule,self.basis,self.method,matrix=self.dmat,multInts=self.multints)
        # multipole moment methods
        self.mtp_methods = []
        # interaction energy methods
        self.eint_methods = []
        #
        self.exchange          = False
        self.has_hexadecapoles = False
        self.cbamm             = False
        self.svd               = False
        self.stat              = False
        self.Print             = False
        self.fit_atoms         = None
        # --- search the tasks
        for task in self.tasks:
            # search the units
            if   'units' in task.lower():
                  self.units = self.search(task)[0]
            # search the unit converter for final result
            if   'result' in task.lower():
                  self.result = self.search(task)[0]
            # exchange 
            if   'exch' in task.lower():
                  self.exchange = self.search(task,bool)
            # search for the interaction energy routines
            if   'eint' in task.lower():
                  self.eint_methods = self.search(task)
            # search for the multipole population routines
            elif 'mtp' in task.lower():
                  self.mtp_methods = self.search(task)
            # search for switching on hexadecapole moments
            if   'hexad' in task.lower():
                  self.has_hexadecapoles = self.search(task,bool)
            # search for switching on CBAMM calculation
            if   'bonds' in task.lower():
                  self.cbamm = self.search(task,bool)
            # search for esp features
            if   'pot' in task.lower():
                  self.pot = self.search(task)[0]
            # search for statistical evaluation of variational space for fitting
            if   'svd' in task.lower():
                  self.svd = self.search(task,bool)
            # search for padding information for ESP fitting
            if   'pad' in task.lower():
                  self.pad = float(self.search(task)[0])
            # search for statistical evaluation
            if   'stat' in task.lower():
                  self.stat = self.search(task,bool)
            # search for printing
            if   'print' in task.lower():
                  self.Print = self.search(task,bool)
            # search for number of points per atom during ESP fit
            if   'mpoints' in task.lower(): 
                  self.mpot = int(self.search(task)[0])
            # search for saving information
            if   'save' in task.lower():
                  self.outname = self.search(task)[0]
            # search for atoms for which ESP has to be performed
            if   'fit_atoms' in task.lower():
                  self.fit_atoms = self.search(task)

        # ---- DO! ----
        # multipole population calculations
        if self.mtp_methods:
           for subtask in self.mtp_methods:
               self.mtp(subtask)
        # coulomb interaction energy calculations
        if self.eint_methods:
           for subtask in self.eint_methods:       
               self.eint(subtask)

    def mtp(self,method):
        """performs multipole population analysis
           routines basing on the procedure (method) and input
           data from PARSER class object."""
           
        if method.lower() == 'camm': 
           print "\n   ------------------------"
           print   "   CAMMS ESTIMATION ROUTINE"
           print   "   ------------------------\n\n\n\n"
           print MULTIP.__printCAMMs__.__doc__
        elif method.lower() == 'mmm':
           print "\n   -----------------------"
           print   "   MMMS ESTIMATION ROUTINE"
           print   "   -----------------------\n\n\n\n"
           print MULTIP.__printMMMs__.__doc__
        elif method.lower() == 'esp':
           print "\n   ------------------------------"
           print   "   ESP CHARGES ESTIMATION ROUTINE"
           print   "   ------------------------------\n\n\n\n"

        if method.lower() != 'esp':
          for i, mol in enumerate(self.M):
            log = "      COMPUTATIONS FOR: --- %s ---" % mol.name
            print "     ","-"*55
            print log
            print "     ","-"*55

            # determine if bond moments should be computed or not
            if not self.cbamm: self.bond_set[i] = None

            # use external multints and density matrix
            try:
               print " READING EXTERNAL DENSITY MATRIX\n"  # EXTERNAL MULTINTS INTERFACE NOT ADDED YET!
               result = MULTIP(mol,self.basis,self.method,
                               matrix=self.dmat_set[i],
                               transition=self.transition_set[i],
                               bonds=self.bond_set[i],
                               hexadecapoles=self.has_hexadecapoles)
            # calculate multints and density matrix using PyQuante
            except IndexError: 
               print " CALCULATION OF DENSITY MATRIX BY DIRECT SCF (OR INDEX ERROR: CHECK IF SCF WAS NOT CHOSEN)"
               result = MULTIP(mol,self.basis,self.method,
                               transition=self.transition_set[i],
                               bonds=self.bond_set[i],
                               hexadecapoles=self.has_hexadecapoles)
               
            if method.lower() == 'camm':
                 result.camms()
                 result.__printCAMMs__()
                 if self.outname is not None:
                    dma = result.get()[0]
                    dma.write(self.outname)
                    self.__dma = dma

            elif method.lower() == 'mmm':
                 result.mmms()
                 #PRINTL(result.P)
                 result.__printMMMs__()
            result.clock.__print__()
            
        else:
           for i, mol in enumerate(self.M):
            log = "      COMPUTATIONS FOR: --- %s ---" % mol.name
            print "     ","-"*55
            print log
            print "     ","-"*55            
            # use external multints and density matrix
            try:
                      print " READING EXTERNAL DENSITY MATRIX\n"  # EXTERNAL MULTINTS INTERFACE NOT ADDED YET!
                      result = ESP(mol,self.basis,self.method,mpot=self.mpot,
                                   pot=self.pot,pad=self.pad,stat=self.stat,SVD=self.svd,
                                   Print=self.Print,
                                   matrix=self.dmat_set[i],
                                   transition=self.transition_set[i],
                                   fit_atoms=self.fit_atoms)
                      if self.fit_atoms is not None:
                         dma = DMA(nfrag=len(self.fit_atoms))
                         dma.set_structure(pos=mol.get_pos()[array(self.fit_atoms,int)-1],equal=True)
                      else:
                         dma = DMA(nfrag=len(mol))
                         dma.set_structure(pos=mol.get_pos()[:]                    ,equal=True)
                      dma.set_moments(charges=result.charges)
                      if self.outname is not None:
                         dma.write(self.outname)
                      self.__dma = dma
                      print result
            except IndexError: 
                      
                      print " CALCULATION OF DENSITY MATRIX BY DIRECT SCF (OR INDEX ERROR: CHECK IF SCF WAS NOT CHOSEN)"
                      result = ESP(mol,self.basis,self.method,mpot=self.mpot,
                                   pot=self.pot,pad=self.pad,stat=self.stat,SVD=self.svd,
                                   Print=self.Print,
                                   transition=self.transition_set[i],
                                   fit_atoms=self.fit_atoms)
                      if self.fit_atoms is not None:
                         dma = DMA(nfrag=len(self.fit_atoms))
                         dma.set_structure(pos=mol.get_pos()[array(self.fit_atoms,int)-1],equal=True)
                      else:
                         dma = DMA(nfrag=len(mol))
                         dma.set_structure(pos=mol.get_pos()[:]                    ,equal=True)
                      dma.set_moments(charges=result.charges)
                      if self.outname is not None:
                         dma.write(self.outname)
                      self.__dma = dma
                      print result
                                   
    def eint(self,method):
        """performs electrostatic interaction energy 
           estimation routines basing on the procedure (method)
           and input data from PARSER class object."""
           
        # for CAMM and ESP interaction energy
        if method.lower() == 'camm': 
           print "\n   -------------------------------------------------------------"
           print   "   MULTIPOLE PART OF ELECTROSTATIC INTERACTION ENERGY FROM CAMMS"
           print   "   -------------------------------------------------------------\n\n\n\n"
        elif method.lower() == 'esp':
           print "\n   -----------------------------------------"
           print   "   ELECTROSTATIC INTERACTION ENERGY FROM ESP"
           print   "   -----------------------------------------\n\n\n\n"

        if method.lower()=='esp' or method.lower()=='camm':
           results = []
           for i,mol in enumerate(self.M):
              text = "      COMPUTATIONS FOR: --- %s ---" % mol.name
              print "     ","-"*55
              print text
              print "     ","-"*55
              if method.lower() == 'camm':
                   # use external multints and density matrix
                   try:
                      result = MULTIP(mol,self.basis,self.method,
                               matrix=self.dmat_set[i],
                               transition=self.transition_set[i],
                               bonds=self.bond_set[i],
                               hexadecapoles=self.has_hexadecapoles)
                      #print "holahola!",self.units
                   # calculate multints and density matrix using PyQuante
                   except IndexError: 
                      result = MULTIP(mol,self.basis,self.method,
                               transition=self.transition_set[i],
                               bonds=self.bond_set[i],
                               hexadecapoles=self.has_hexadecapoles)
                   result.camms()
                   result.makeTracelessCAMMs()
                   result.__printCAMMs__()
              elif method.lower() == 'esp':
                   # use external multints and density matrix
                   try:
                      result = ESP(mol,self.basis,self.method,mpot=self.mpot,
                                   pot=self.pot,pad=self.pad,stat=self.stat,SVD=self.svd,
                                   Print=self.Print,
                                   matrix=self.dmat_set[i],
                                   transition=self.transition_set[i],
                                   fit_atoms=self.fit_atoms) 
                   except IndexError: 
                      result = ESP(mol,self.basis,self.method,mpot=self.mpot,
                                   pot=self.pot,pad=self.pad,stat=self.stat,SVD=self.svd,
                                   Print=self.Print,
                                   transition=self.transition_set[i],
                                   fit_atoms=self.fit_atoms) 
              results.append(result)
              #result.clock.__print__()

           if  method.lower() == 'camm':
               etext = 'DE(EL,MTP,CAMM)'
               Eint = Eelcamm(results[0],results[1])
               print Eelcamm.log
           elif method.lower() == 'esp':
               etext = 'DE(EL,ESP)'
               Eint = EelESP(results[0],results[1])      

           print 
           print " %s = %12.8f [CM-1]" % (etext,Eint*UNITS.HartreeToCmRec)
           print
           results[1].clock.actualize('%s interaction energy' % etext)
           results[1].clock.__print__()
         
        # for EEL EDS interaction energy
        if method.lower() == 'eeleds':
           print "\n   -----------------------------------------"
           print   "   ELECTROSTATIC INTERACTION ENERGY FROM EDS"
           print   "   -----------------------------------------\n\n\n\n"
           try:
              result = EELEDS(self.M[0],self.M[1],
                         self.basis,self.method,
                         matrixa=self.dmat_set[0],
                         matrixb=self.dmat_set[1],
                         transition=self.transition_set[0],
                         exchange=self.exchange)
           except IndexError: 
              result = EELEDS(self.M[0],self.M[1],
                         self.basis,self.method,
                         transition=self.transition_set[0],
                         exchange=self.exchange)
           print result

    def search(self,task,dtype=str):
        """Withdraws the subtask. Format in the input file:
TASK=subtask1,subtask2[,...]"""
        t = task.split('=')[1].split(',')
        if   dtype==str:                   return t
        elif dtype==bool:
             if   t[0].lower() == 'true' : return True
             elif t[0].lower() == 'false': return False
             elif t[0].lower() == '1'    : return True
             elif t[0].lower() == '0'    : return False
             else: raise ValueError, "Incorrect value for task %s" % task


