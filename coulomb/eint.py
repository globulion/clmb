# --------------------------------------------------- #
#        COULOMB INTERACTION ENERGY EVALUATORS        #
# --------------------------------------------------- #

from libbbg.units import *
import numpy

__all__=['Eeldc','EelESP','Eelcamm']

def Eeldc(dc1,dc2):
    """calculates the coulombic part of interaction 
       energy from two DC (Density Cubes) with data
       memorial in memory (discalloc for both cube
       file objects has to be set off)""" 
       
    Eint=0
    # origins of cubes
    o1 = dc1.origin; o2 = dc2.origin
    # cube dimensions 
    nxa, nya, nza = dc1.dimension
    nxb, nyb, nzb = dc2.dimension
    # cube spacings
    vxa, vya, vza = dc1.spacings
    vxb, vyb, vzb = dc2.spacings
    #
    for i,a in enumerate(dc1.data):
        for j,b in enumerate(dc2.data):
            # point-at-once indexing
            ixa = i/(nya*nza)
            iya = (i-ixa*nya*nza)/nza
            iza = i-ixa*nya*nza-iya*nza
            #
            ixb = i/(nyb*nzb)
            iyb = (i-ixb*nyb*nzb)/nzb
            izb = i-ixb*nyb*nzb-iyb*nzb
            # point-at-once position
            ri = o1 + vxa*ixa + vya*iya + vza*iza 
            rj = o2 + vxb*ixb + vyb*iyb + vzb*izb 
            # 
            rij = numpy.sqrt(numpy.sum((ri-rj)**2)) 
            Eint+= a*b/rij
            #
    return Eint        

def EelESP(esp1,esp2):
    """calculates E(EL) from two charge distributions 
esp1 and esp2 obtained from fitting these charges to
the molecular electrostatic potential. Energy is given
in a.u."""

    q1,q2 = esp1.charges, esp2.charges
    r1,r2 = esp1.RArray,  esp2.RArray
    N1,N2 = len(q1),len(q2)
    Eint=0
    for i in range(N1):
        for j in range(N2):
            Eint+=q1[i]*q2[j]/numpy.sqrt(numpy.sum((r1[i]-r2[j])**2))
    return Eint

def Eelcamm(camm1,camm2,converter=UNITS.HartreePerHbarToCmRec):
    """calculates E(EL)MTP from two CAMM distributions.
camm1 and camm2 are the objects of the class MULT after 
performing the 'camms()' method on them. Calculations are
in atomic units and a respective interaction energy is in 
a.u. as well. Remember to convert CAMMs to traceless form!
('makeTracelessCAMMs()' method)"""

    #camm1.makeTracelessCAMMs()
    #camm2.makeTracelessCAMMs()
    if camm1.has_hexadecapoles and camm2.has_hexadecapoles:
       hexadecapoles = True
       Ra,qa,Da,Qa,Oa,Ha = camm1.ReturnCAMMs()
       Rb,qb,Db,Qb,Ob,Hb = camm2.ReturnCAMMs()
    else:
       hexadecapoles = False
       Ra,qa,Da,Qa,Oa    = camm1.ReturnCAMMs()
       Rb,qb,Db,Qb,Ob    = camm2.ReturnCAMMs()
    #
    qq = 0
    qD = 0 ; Dq = 0
    qQ = 0 ; Qq = 0
    DQ = 0 ; QD = 0
    QQ = 0 
    qO = 0 ; Oq = 0
    DD = 0 
    DO = 0 ; OD = 0
    QO = 0 ; OQ = 0
    OO = 0 ;
    qH = 0 ; Hq = 0
    #
    Tensordot = numpy.tensordot
    for i in xrange(len(Ra)):
         for j in xrange(len(Rb)):
             R    = Rb[j]-Ra[i]
             Rab=numpy.sqrt(numpy.sum(R**2,axis=0))
             qq  +=   qa[i]*qb[j]/Rab                                                               # qa - qb  | R1
             qD  +=  -qa[i]*Tensordot(Db[j],R,(0,0))/Rab**3                                         # qa - Db  | R2
             Dq  +=  +qb[j]*Tensordot(Da[i],R,(0,0))/Rab**3                                         # qb - Da  | R2
             DD  +=-3*Tensordot(Da[i],R,(0,0))*Tensordot(Db[j],R,(0,0))/Rab**5                      # Da - Db  | R3
             DD  +=   Tensordot(Da[i],Db[j],(0,0))/Rab**3                                           # Da - Db  | R3
             qQ  +=   qa[i]*Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                      # qa - Qb  | R3
             Qq  +=   qb[j]*Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                      # qb - Qa  | R3
             DQ  +=-2*Tensordot(Da[i],Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                        # Da - Qb  | R4
             QD  += 2*Tensordot(Db[j],Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                        # Db - Qa  | R4
             DQ  += 5*Tensordot(Da[i],R,(0,0))*Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))/Rab**7   # Da - Qb  | R4
             QD  +=-5*Tensordot(Db[j],R,(0,0))*Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))/Rab**7   # Db - Qa  | R4
             qO  +=  -qa[i]*Tensordot(R,Tensordot(R,Tensordot(Ob[j],R,(0,0)),(0,0)),(0,0))/Rab**7   # qa - Ob  | R4
             Oq  +=   qb[j]*Tensordot(R,Tensordot(R,Tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7   # qb - Oa  | R4
             QQ  += (35.)/(3.)* (Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0)) *
                                 Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))  ) / Rab**9            # Qa - Qb  | R5
             OD  +=-7*(Tensordot(Db[j],R,(0,0)) *
                       Tensordot(R,Tensordot(R,Tensordot(Oa[i],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Db - Oa  | R5
             DO  +=-7*(Tensordot(Da[i],R,(0,0)) *
                       Tensordot(R,Tensordot(R,Tensordot(Ob[j],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Da - Ob  | R5
             QQ  +=-(20.)/(3.) * Tensordot(Tensordot(R,Qa[i],(0,0)),
                                           Tensordot(R,Qb[j],(0,0)),(0,0)) / Rab**7                 # Qa - Qb  | R5
             QQ  +=(2.)/(3.)  * Tensordot(Qa[i],Qb[j])  / Rab**5                                    # Qa - Qb  | R5
             OD  +=3 * Tensordot(R,Tensordot(R,Tensordot(Oa[i],Db[j],(0,0)),(0,0)),(0,0)) / Rab**7  # Db - Oa  | R5
             DO  +=3 * Tensordot(R,Tensordot(R,Tensordot(Ob[j],Da[i],(0,0)),(0,0)),(0,0)) / Rab**7  # Da - Ob  | R5
             if hexadecapoles:
                Hq  += qb[j] * Tensordot(R,Tensordot(R,Tensordot(R,Tensordot(R,Ha[i],                               
                                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                              # Ha - qb  | R5
                qH  += qa[i] * Tensordot(R,Tensordot(R,Tensordot(R,Tensordot(R,Hb[j],
                                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                              # Hb - qj  | R5
             ### these are implemented already !
             OQ  += 2* Tensordot(Tensordot(Oa[i],Qb[j],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qb - Oa  | R6
             QO  +=-2* Tensordot(Tensordot(Ob[j],Qa[i],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qa - Ob  | R6
             OQ  +=-14*Tensordot(Tensordot(R,Tensordot(Oa[i],R,(1,0)),(0,0)) ,                      # Qb - Oa  | R6
                                 Tensordot(R,Qb[j],(0,0)) ,(0,0)) / Rab**9                          
             QO  += 14*Tensordot(Tensordot(R,Tensordot(Ob[j],R,(1,0)),(0,0)) ,                      # Qa - Ob  | R6
                                 Tensordot(R,Qa[i],(0,0)) ,(0,0)) / Rab**9
             OQ  +=( 21*Tensordot(Tensordot(R,Tensordot(Oa[i],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
                      * Tensordot(R,Tensordot(Qb[j],R,(0,0)),(0,0))) / Rab**11
             QO  +=(-21*Tensordot(Tensordot(R,Tensordot(Ob[j],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
                      * Tensordot(R,Tensordot(Qa[i],R,(0,0)),(0,0))) / Rab**11   
             OO  +=(2.)/(5.)*Tensordot(Oa[i],Ob[j],((0,1,2),(0,1,2))) / Rab**7                      # Ob - Oa  | R7
             OO  +=(-42./5.)*Tensordot(Tensordot(R,Oa[i],(0,0)),
                                       Tensordot(R,Ob[j],(0,0)),
                                       ((0,1),(0,1))) / Rab**9                                      # Ob - Oa  | R7
             OO  +=(189.)/(5.)*Tensordot(
                                         Tensordot(Tensordot(R,Oa[i],(0,0)),R,(0,0)),
                                         Tensordot(Tensordot(R,Ob[j],(0,0)),R,(0,0)),
                                         (0,0)) /Rab**11                                            # Oa - Ob  | R7
             OO  +=-(231./5.)*(Tensordot(Tensordot(Tensordot(R,Oa[i],(0,0)),R,(0,0)),R,(0,0)) *
                               Tensordot(Tensordot(Tensordot(R,Ob[j],(0,0)),R,(0,0)),R,(0,0)) ) /\
                               Rab**13                                                              # Oa - Ob  | R7
            
             Eint = qq + qD + Dq + qQ + Qq + qO + Oq + DD + DQ + QD + DO + OD + QQ + QO + OQ + OO + qH + Hq
             
             ### save the partitioning for current usage
             Eelcamm.qq = qq;Eelcamm.qD = qD;Eelcamm.qQ = qQ;Eelcamm.qO = qO;Eelcamm.QO = QO;Eelcamm.qH = qH
             Eelcamm.DD = DD;Eelcamm.DQ = DQ;Eelcamm.DO = DO;Eelcamm.QQ = QQ;Eelcamm.OO = OO;Eelcamm.Hq = Hq
    #
    Eelcamm.A = (            qq )
    Eelcamm.B = (Eelcamm.A + qD + Dq )
    Eelcamm.C = (Eelcamm.B + DD + qQ + Qq )
    Eelcamm.D = (Eelcamm.C + qO + Oq + DQ + QD)
    Eelcamm.E = (Eelcamm.D + DO + OD + QQ + qH + Hq)
    Eelcamm.F = (Eelcamm.A + Dq + qD + DD )
    Eelcamm.G = (Eelcamm.F + Qq + qQ + QD + DQ + QQ)
    Eelcamm.H = (Eelcamm.G + Oq + qO + OD + DO + OQ + QO + OO)
    #
    Eelcamm.A *= converter
    Eelcamm.B *= converter
    Eelcamm.C *= converter
    Eelcamm.D *= converter
    Eelcamm.E *= converter
    #
    log = "\n" 
    log+= " --------------------------------:--------------------------\n"
    log+= " INTERACTION ENERGY TERMS [cm-1] : PARTITIONING TERMS [cm-1]\n"
    log+= " --------------------------------:--------------------------\n"
    log+= "%6s %20.2f      :\n" % ("Total".rjust(6),Eint*converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "%7s %19.2f      :  1        %10.2f\n" % ("q-q".rjust(6), qq    *converter,Eelcamm.A)
    log+= "%7s %19.2f      :  1+2      %10.2f\n" % ("q-D".rjust(6),(qD+Dq)*converter,Eelcamm.B)
    log+= "%7s %19.2f      :  1+2+3    %10.2f\n" % ("q-Q".rjust(6),(qQ+Qq)*converter,Eelcamm.C)
    log+= "%7s %19.2f      :  1+2+3+4  %10.2f\n" % ("q-O".rjust(6),(qO+Oq)*converter,Eelcamm.D)
    log+= "%7s %19.2f      :  1+2+3+4+5%10.2f\n" % ("D-D".rjust(6), DD    *converter,Eelcamm.E)
    log+= "%7s %19.2f      :\n"                  % ("D-Q".rjust(6),(DQ+QD)*converter)
    log+= "%7s %19.2f      :\n"                  % ("D-O".rjust(6),(DO+OD)*converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-Q".rjust(6), QQ    *converter)
    log+= "%7s %19.2f      :\n"                  % ("Q-O".rjust(6),(QO+OQ)*converter)
    log+= "%7s %19.2f      :\n"                  % ("O-O".rjust(6), OO    *converter)
    log+= "%7s %19.2f      :\n"                  % ("q-H".rjust(6),(qH+Hq)*converter)
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "\n"
    Eelcamm.log = log

    return Eint


# --------------- Interaction energy ---------------------- # 
