#!/usr/bin/python
# *-* coding:utf-8 *-*
"""
 COULOMB Introductory information
"""
# -------------------------------------
from setuptools import setup, Extension
# -------------------------------------
import sys
# --------

def check_python_version():
    python_version = sys.version_info
    print(" Python version is ", python_version[:3])
    if not python_version[:2] in [(2, 7), (2, 6)]: 
       print(" VersionError: Coulomb requires Python 2.6 or Python 2.7!")
       sys.exit(-1)
    return

def Main(argv):
    check_python_version()
    classifiers = ['Development Status :: 1 - Planning',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: Freeware',
                   'Operating System :: POSIX',
                   'Operating System :: MacOS :: MacOS X',
                   'Operating System :: Unix',
                   'Programming Language :: Python',
                   'Programming Language :: C',
                   'Programming Language :: Fortran',
                   'Topic :: Scientific/Engineering :: Bio-Informatics',
                   'Topic :: Scientific/Engineering :: Chemistry',
                   'Topic :: Scientific/Engineering :: Physics',
                   'Topic :: Software Development :: Libraries :: Python Modules',
                   ]
    extensions  = []
    print(" Installing has started!")
    setup(name                         =  'Coulomb_t'  , 
          version                      =  '0.0.0'     ,
          description                  =  'The program for electrostatic quantum chemistry computations' ,
          long_description             =   open("README.md").read(),
          author                       =  'Bartosz Błasiak'      ,
          author_email                 =  'blasiak.bartosz@gmail.com',
          contributors                 = ['Robert Władysław Góra',],
          contributors_emails          = ['robert.gora@pwr.edu.pl',],
          url                          =  '',
          license                      =   open("LICENSE.md").read(),
          requires                     = ['numpy (>=1.0.3)', 'PyQuante',],
          packages                     = ['Coulomb_t',
                                          'Coulomb_t.core',
                                          'Coulomb_t.qm',
                                          'Coulomb_t.data',],
          package_dir                  = {'Coulomb_t': 'Coulomb_t'},
          classifiers                  =   classifiers,
          ext_modules                  =   extensions,
          install_requires             = ['numpy>=1.4.3',
                                          'PyQuante',],
          test_suite                   =  '',
          tests_require                = [],
          )
          

if __name__ == "__main__": Main(sys.argv[1:])


