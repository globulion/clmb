﻿# --------------------------------------------- #
#               WORKING UTILITIES               #
# --------------------------------------------- #

from units  import *
from parser import PARSER
from multip import * 
from esp    import *
from eeleds import *
from eint   import *  

__all__=['DO']

class DO(PARSER):
    """contains main routine for COULOMB.py"""
    
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
        # --- search the tasks
        for task in self.tasks:
            # search the units
            if   'units' in task.lower():
                  self.units = self.search(task)[0]
            # search the unit converter for final result
            if   'result' in task.lower():
                  self.result = self.search(task)[0]
            # transition 
            if   'trans' in task.lower():
                  self.transition = True
            # state (if transition)
            if   'trst' in task.lower():
                  self.state = self.search(task)[0]
            # exchange 
            if   'exch' in task.lower():
                  self.exchange = True  
            # search for the interaction energy routines
            if   'eint' in task.lower():
                  self.eint_methods = self.search(task)
            # search for the multipole population routines
            elif 'mtp' in task.lower():
                  self.mtp_methods = self.search(task)
            # search for esp features
            if   'pot' in task.lower():
                  self.pot = self.search(task)[0]
            if   'svd' in task.lower():
                  self.svd = bool(self.search(task)[0])
            if   'pad' in task.lower():
                  self.pad = float(self.search(task)[0])
            if   'stat' in task.lower():
                  self.stat = bool(self.search(task)[0])
            if   'print' in task.lower():
                  self.Print = bool(self.search(task)[0])
            if   'mpoints' in task.lower(): 
                  self.mpot = int(self.search(task)[0])
                  
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
            # use external multints and density matrix
            try:
               result = MULTIP(mol,self.basis,self.method,
                               matrix=self.dmat_set[i],
                               transition=self.transition)
            # calculate multints and density matrix using PyQuante
            except IndexError: 
               result = MULTIP(mol,self.basis,self.method,
                               transition=self.transition)
               
            if method.lower() == 'camm':
                 result.camms()
                 result.__printCAMMs__()

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
                      result = ESP(mol,self.basis,self.method,mpot=self.mpot,
                                   pot=self.pot,pad=self.pad,stat=self.stat,SVD=self.svd,
                                   Print=self.Print,
                                   matrix=self.dmat_set[i],
                                   transition=self.transition) 
            except IndexError: 
                      result = ESP(mol,self.basis,self.method,mpot=self.mpot,
                                   pot=self.pot,pad=self.pad,stat=self.stat,SVD=self.svd,
                                   Print=self.Print,
                                   transition=self.transition) 
                                   
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
                               transition=self.transition)
                      #print "holahola!",self.units
                   # calculate multints and density matrix using PyQuante
                   except IndexError: 
                      result = MULTIP(mol,self.basis,self.method,
                               transition=self.transition)
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
                                   transition=self.transition) 
                   except IndexError: 
                      result = ESP(mol,self.basis,self.method,mpot=self.mpot,
                                   pot=self.pot,pad=self.pad,stat=self.stat,SVD=self.svd,
                                   Print=self.Print,
                                   transition=self.transition) 
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
                         transition=self.transition,
                         exchange=self.exchange)
           except IndexError: 
              result = EELEDS(self.M[0],self.M[1],
                         self.basis,self.method,
                         transition=self.transition,
                         exchange=self.exchange)
           print result

    def search(self,task):
        """withdraws the subtask. Format in the input file:
               TASK=subtask1,subtask2[,...]"""
        return task.split('=')[1].split(',')


