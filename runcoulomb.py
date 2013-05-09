#!/usr/bin/python
# -*- coding: utf-8 -*-

#     Copyright © 2012, Bartosz Błasiak (globula@o2.pl)
#                       Robert W. Góra  (robert.gora@pwr.wroc.pl)
#  
#     This program is assigned for free use. You can modify the code
#     as you wish. If any bugs found please send them to me via mail.
#
#     The script was prepared for the purpose of computational tasks
#     specified in grant no. 

# ----------------------------------------------------------------------------------------------

from coulomb_head import *

# ----------------- #
#    MAIN ROUTINE   #
# ----------------- #

def Main(argv):
    """ Coulomb.py main routine """
    
    # THE DEFAULT PARAMETERS
    input = ''
    
    # OPTIONS AND FUNCTIONALITY MAP 
    try:
        opts, args = getopt.getopt(argv, "hvf:"       ,
                                        ["help"       ,
                                         "version"    ])
    except getopt.GetoptError, error: Error()
    if not argv: Usage()
    
    for opt, arg in opts:
        # --- help
        if opt in ("-h", "--help"    ):
           Usage()
        # --- version info
        if opt in ("-v", "--version" ):
           Version()
        # --- read input file
        
    if args:
       input = args[0]
       print __doc__
       DO(input)
       
    # print the input file
    if input:
       print "\n\n ==================================="
       print     " ECHO OF INPUT FILE:\n"
       for line in open(input).readlines():
           print ' >> INPUT CARD: ', line,
       print


# ----------------- #
#        DO !       #
# ----------------- #

# ---------------------------------------
if __name__ == '__main__': Main(argv[1:])
# ---------------------------------------
#
#                  Globulion@je@bulion@
# 
#cube_test()
