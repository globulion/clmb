# ------------------------------------- #
#               T E S T S               #
# ------------------------------------- #

from coulomb_head import *

def cube_test():
    """testy cubów"""
    dc1 = CUBFLE('da.cube')
    dc2 = CUBFLE('cat.cube')
    print Eeldc(dc1,dc2)
    pass
def i():
 """testy"""
 o_cyjano = Molecule('o-CYJANOFENOL',

 [(6,       (-6.51624,       3.68194,       0.06330)),
  (6,       (-5.46325,       2.76963,       0.16100)),
  (6,       (-4.13500,       3.21623,       0.12648)),
  (6,       (-3.87506,       4.58499,      -0.00682)),
  (6,       (-6.24765,       5.04346,      -0.06922)),
  (6,       (-4.92802,       5.49431,      -0.10422)),
  (6,       (-3.07564,       2.25975,       0.22858)),
  (1,       (-7.06570,       5.75597,      -0.14549)),
  (1,       (-7.54574,       3.33154,       0.09048)),
  (1,       (-5.68718,       1.70937,       0.26409)),
  (8,       (-2.61769,       5.12452,      -0.05103)),
  (1,       (-4.71908,       6.55621,      -0.20759)),
  (1,       (-1.93672,       4.43303,       0.02217)),
  (7,       (-2.18806,       1.51648,       0.30839))],
  units='Angstrom')

 o_cyjano = Molecule('o-CYJANOFENOL',

 [(6,       (-6.71440,       3.74633,       0.05498 )),
  (6,       (-5.71653,       2.77876,       0.18804 )),
  (6,       (-4.36191,       3.14689,       0.18977 )),
  (6,       (-4.02000,       4.49603,       0.05069 )),
  (6,       (-6.36563,       5.08840,      -0.07428 )),
  (6,       (-5.02132,       5.46131,      -0.07534 )),
  (6,       (-3.37848,       2.11942,       0.34741 )),
  (1,       (-7.13826,       5.84680,      -0.17442 )),
  (1,       (-7.76174,       3.45376,       0.05493 )),
  (1,       (-6.00344,       1.73366,       0.29318 )),
  (8,       (-2.74196,       4.98204,       0.03085 )),
  (1,       (-4.74938,       6.50946,      -0.17632 )),
  (1,       (-2.05896,       4.27527,      -0.01872 )),
  (7,       (-2.56924,       1.29389,       0.46977 ))],
  units='Angstrom')


 ch3oh = Molecule('CH3OH',
[(8,         (-0.66105,       5.33623,      -0.19272)),
 (1,         (-0.74525,       4.40889,       0.11994)),
 (6,         ( 0.50624,       5.36701,      -0.99930)),
 (1,         ( 0.63392,       6.37510,      -1.39991)),
 (1,         ( 0.39591,       4.65927,      -1.82462)),
 (1,         ( 1.37516,       5.10344,      -0.39155))],
units='Angstrom')

 h2o = Molecule('H2O',
                    [(8,  ( 1.00000000,     0.00000000,     0.04851804)),
                     (1,  ( 1.75300223,     0.20000000,    -0.51923377)),
                     (1,  (-0.75300223,     0.00000000,    -0.51923377))],
                    units='Angstrom')

 h2 = Molecule('H2' ,
                    [(1,  ( 1.75300223,     3.20000000,    -1.51923377)),
                     (1,  (-0.75300223,     3.00000000,    -1.51923377))],
                    units='Angstrom')

 nh3 = Molecule('NH3',
[(7,       (  1.74255,       1.66277,       0.64019)),
 (1,       (  2.65514,       1.28877,       0.36292)),
 (1,       (  1.86554,       2.67852,       0.62003)),
 (1,       (  1.66089,       1.42925,       1.63284))],units='Angstrom')
 h2o=Molecule('H2O',
[(8,       (  0.00560,       0.53892,      -1.24338)),
 (1,       (  0.67353,       0.91132,      -0.62397)),
 (1,       ( -0.82070,       0.80348,      -0.81432))],units='Angstrom')


 #camm1 = MULTIP(h2o,argv[1],'HF')
 #camm1.camms()
 #camm1.mmms()
 #camm1.__printCAMMs__()
 #camm1.__printMMMs__()
 #camm1.clock.__print__() 
# MMM=int(argv[1])
 #g = RUN(h2o,'sto3g','HF')
 #h = MULTIP(h2o,'sto-3g','HF')
 #h.camms()
 #h.makeTracelessCAMMs()
 #R=array([2.33,3.2334,-4.34])
 #V1 = Vr_wfn(g,R)
 #V2 = Vr_camms(h,R)
 #print V2
 #print V1
 basis='6-31G**';method='HF'
 basis='6-31G'
 basis='sto3g'
 mol1 = h2o
 mol2 = h2
# r=DFI(h2,basis,method,h2o,basis,method,conv=1e-10)
 r = EELEDS(mol1,mol2,basis,'HF',Exchange=False)
 print 
 #
 t1 = MULTIP(mol1,basis,'HF')
 t1.camms()
 t1.makeTracelessCAMMs() 
 t2 = MULTIP(mol2,basis,'HF')
 t2.camms()
 t2.makeTracelessCAMMs() 
 print Eelcamm(t2,t1)
 print Eelcamm(t1,t2)
 #p1=read_matrix('gaussian',argv[1],'Beta transition density to state',1)
 #p1+=read_matrix('gaussian',argv[1],'Alpha transition density to state',1)
 #p2=read_matrix('gaussian',argv[2],'Beta transition density to state',1)
 #p2+=read_matrix('gaussian',argv[2],'Alpha transition density to state',1)
 #print " Transition density matrices derived. Now 2-el integrals BTF ..."

 e1= Molecule('ethylene1' ,
[(6,                    (-0.00010, -0.70254,  2.50001)),
 (1,                    (-0.91646, -1.25118,  2.60255)),
 (6,                    (-0.0001 ,  0.70255,  2.50000)),
 (1,                    ( 0.91712, -1.25073,  2.60004)),
 (1,                    ( 0.91713,  1.25065,  2.59922)),
 (1,                    (-0.91663,  1.2512 ,  2.60269))],
                    units='Angstrom')
 e2= Molecule('ethylene2' ,
[(6,                   (-0.0001 , -0.70254, -0.50001)),
 (1,                   (-0.91646, -1.25118, -0.60255)),
 (6,                   (-0.0001 ,  0.70255, -0.5    )),
 (1,                   ( 0.91712, -1.25073, -0.60004)),
 (1,                   ( 0.91713,  1.25065, -0.59922)),
 (1,                   (-0.91663,  1.2512 , -0.60269))],
                    units='Angstrom')


 b='6-31G(d)'
 #TDFI(e1,b,p1,e2,b,p2)
 #t1 = MULTIP(e1,b,'HF',matrix=p1,Transition=True)
 #t1.camms()
 #t1.makeTracelessCAMMs() 
 #t2 = MULTIP(e2,b,'HF',matrix=p2,Transition=True)
 #t2.camms()
 #t2.makeTracelessCAMMs() 
 Bohr_To_CmInv = 219474.63
 #print Eelcamm(t2,t1)*Bohr_To_CmInv
 #print Eelcamm(t1,t2)*Bohr_To_CmInv
 #TDFI(e1,b,p1,e2,b,p2)

 #for i in range(38):
 #    for j in range(38):
 #        print p1[i,j]-p2[i,j], 
 #    print
 #  
 #t1 = ESP(mol1,basis,'HF',mpot=5000*len(mol1.atoms),pot='WFN',SVD=False,stat=False)
 #t2 = ESP(mol2,basis,'HF',mpot=5000*len(mol2.atoms),pot='WFN',SVD=False,stat=False)
 #print EelESP(t1,t2)
# v =ESP(h2o,'sto-3g','HF',mpot=MMM,pot='CAMM',stat=True,SVD=True) 
# v.clock.__print__()
# r = open('test.dat','w')
# for i in range(MMM):
#     #print >> r, v.VArray[i][0], v.VArray[i][1], v.VArray[i][3]
#     y=(v.VFArray[i][3]-v.VArray[i][3])/v.VArray[i][3]*100
#     print >> r, "%10.3f %10.3f %10.3f %12.8f %12.8f %10.1f%%"%(v.VArray[i][0], v.VArray[i][1], v.VArray[i][2], v.VArray[i][3], v.VFArray[i][3], y)
#print 
#i()
