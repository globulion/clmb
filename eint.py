# --------------------------------------------------- #
#        COULOMB INTERACTION ENERGY EVALUATORS        #
# --------------------------------------------------- #

from units        import *
from numpy import *

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
            rij = sqrt(sum((ri-rj)**2)) 
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
            Eint+=q1[i]*q2[j]/sqrt(sum((r1[i]-r2[j])**2))
    return Eint

def Eelcamm(camm1,camm2):
    """calculates E(EL)MTP from two CAMM distributions.
camm1 and camm2 are the objects of the class MULT after 
performing the 'camms()' method on them. Calculations are
in atomic units and a respective interaction energy is in 
a.u. as well. Remember to convert CAMMs to traceless form!
('makeTracelessCAMMs()' method)"""

    #camm1.makeTracelessCAMMs()
    #camm2.makeTracelessCAMMs()
    Ra,qa,Da,Qa,Oa = camm1.ReturnCAMMs() # hexadecapole integrals not implemented yet
    Rb,qb,Db,Qb,Ob = camm2.ReturnCAMMs()
    #
    converter=UNITS.HartreePerHbarToCmRec
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
    for i in range(len(Ra)):
         for j in range(len(Rb)):
             R    = Rb[j]-Ra[i]
             Rab=sqrt(sum(R**2,axis=0))
             qq  +=   qa[i]*qb[j]/Rab                                                               # qa - qb  | R1
             #if not hash:
             qD  +=  -qa[i]*tensordot(Db[j],R,(0,0))/Rab**3                                         # qa - Db  | R2
             Dq  +=  +qb[j]*tensordot(Da[i],R,(0,0))/Rab**3                                         # qb - Da  | R2
             DD  +=-3*tensordot(Da[i],R,(0,0))*tensordot(Db[j],R,(0,0))/Rab**5                      # Da - Db  | R3
             DD  +=   tensordot(Da[i],Db[j],(0,0))/Rab**3                                           # Da - Db  | R3
             qQ  +=   qa[i]*tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                      # qa - Qb  | R3
             Qq  +=   qb[j]*tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                      # qb - Qa  | R3
             DQ  +=-2*tensordot(Da[i],tensordot(Qb[j],R,(0,0)),(0,0))/Rab**5                        # Da - Qb  | R4
             QD  += 2*tensordot(Db[j],tensordot(Qa[i],R,(0,0)),(0,0))/Rab**5                        # Db - Qa  | R4
             DQ  += 5*tensordot(Da[i],R,(0,0))*tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))/Rab**7   # Da - Qb  | R4
             QD  +=-5*tensordot(Db[j],R,(0,0))*tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))/Rab**7   # Db - Qa  | R4
             qO  +=  -qa[i]*tensordot(R,tensordot(R,tensordot(Ob[j],R,(0,0)),(0,0)),(0,0))/Rab**7   # qa - Ob  | R4
             Oq  +=   qb[j]*tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7   # qb - Oa  | R4
             QQ  += (35.)/(3.)* (tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0)) *
                                 tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))  ) / Rab**9            # Qa - Qb  | R5
             OD  +=-7*(tensordot(Db[j],R,(0,0)) *
                       tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Db - Oa  | R5
             DO  +=-7*(tensordot(Da[i],R,(0,0)) *
                       tensordot(R,tensordot(R,tensordot(Ob[j],R,(0,0)),(0,0)),(0,0)) ) / Rab**9    # Da - Ob  | R5
             QQ  +=-(20.)/(3.) * tensordot(tensordot(R,Qa[i],(0,0)),
                                           tensordot(R,Qb[j],(0,0)),(0,0)) / Rab**7                 # Qa - Qb  | R5
             QQ  +=(2.)/(3.)  * tensordot(Qa[i],Qb[j])  / Rab**5                                    # Qa - Qb  | R5
             OD  +=3 * tensordot(R,tensordot(R,tensordot(Oa[i],Db[j],(0,0)),(0,0)),(0,0)) / Rab**7  # Db - Oa  | R5
             DO  +=3 * tensordot(R,tensordot(R,tensordot(Ob[j],Da[i],(0,0)),(0,0)),(0,0)) / Rab**7  # Da - Ob  | R5
             ### The remaining terms with hexadecapoles are not implemented yet
             #Eint+= qb[j] * tensordot(R,tensordot(R,tensordot(R,tensordot(R,Ha[i],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Ha - qb  | R5
             #Eint+= qa[i] * tensordot(R,tensordot(R,tensordot(R,tensordot(R,Hb[j],
             #                (0,0)),(0,0)),(0,0)),(0,0))   / Rab**9                                 # Hb - qj  | R5
             ### these are implemented already !
             #OQ  += 2* tensordot(tensordot(Oa[i],Qb[j],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qb - Oa  | R6
             #QO  +=-2* tensordot(tensordot(Ob[j],Qa[i],((0,1),(0,1))),R,(0,0)) / Rab**7             # Qa - Ob  | R6
             #OQ  +=-14*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)) ,                      # Qb - Oa  | R6
             #                    tensordot(R,Qb[j],(0,0)) ,(0,0)) / Rab**9                          
             #QO  += 14*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)) ,                      # Qa - Ob  | R6
             #                    tensordot(R,Qa[i],(0,0)) ,(0,0)) / Rab**9
             #OQ  +=( 21*tensordot(tensordot(R,tensordot(Oa[i],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * tensordot(R,tensordot(Qb[j],R,(0,0)),(0,0))) / Rab**11
             #QO  +=(-21*tensordot(tensordot(R,tensordot(Ob[j],R,(1,0)),(0,0)),R,(0,0))              # Qb - Oa  | R6
             #         * tensordot(R,tensordot(Qa[i],R,(0,0)),(0,0))) / Rab**11   
             #OO  +=(2.)/(5.)*tensordot(Oa[i],Ob[j],((0,1,2),(0,1,2))) / Rab**7                      # Ob - Oa  | R7
             #OO  +=(-42./5.)*tensordot(tensordot(R,Oa[i],(0,0)),
             #                          tensordot(R,Ob[j],(0,0)),
             #                          ((0,1),(0,1))) / Rab**9                                      # Ob - Oa  | R7
             #OO  +=(189.)/(5.)*tensordot(
             #                            tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),
             #                            tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),
             #                            (0,0)) /Rab**11
             #OO  +=-(231./5.)*(tensordot(tensordot(tensordot(R,Oa[i],(0,0)),R,(0,0)),R,(0,0)) *
             #                  tensordot(tensordot(tensordot(R,Ob[j],(0,0)),R,(0,0)),R,(0,0)) ) /\
             #                  Rab**13
             
             Eint = qq + qD + Dq + qQ + Qq + qO + Oq + DD + DQ + QD + DO + OD + QQ + QO + OQ + OO
             
             ### save the partitioning for current usage
             Eelcamm.qq = qq;Eelcamm.qD = qD;Eelcamm.qQ = qQ;Eelcamm.qO = qO;Eelcamm.QO = QO;
             Eelcamm.DD = DD;Eelcamm.DQ = DQ;Eelcamm.DO = DO;Eelcamm.QQ = QQ;Eelcamm.OO = OO;
    #
    Eelcamm.A = (         qq )
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
    log+= " "+"-"*32+":"+"-"*26+"\n"
    log+= "\n"
    Eelcamm.log = log

    return Eint


# --------------- Interaction energy ---------------------- # 
