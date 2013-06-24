# ------------------------------------------------ #
#        ELECTROSTATIC POTENTIAL EVALUATORS        #
# ------------------------------------------------ #

__all__=['Vr_wfn','Vr_wfn_1','Vr_camms']

from numpy import *

def Vr_wfn(run,Rb):
    """calculates electrostatic potential in point Rb
directly from wave function. The first argument is an 
object of the 'RUN' class. Rb is in AU."""

    a = run.molecule.atoms
    b = run.bfs
    P = run.P
    Rb= tuple(Rb)
    V = 0
    # nuclear contribution
    for q in range(run.N_at):
        V+= a[q].atno / sqrt( sum( (array(a[q].pos())-Rb)**2 ) )
    # electronic contribution 
    for i in range(run.K):
        for j in range(run.K):
            V += P[i,j] * b[i].nuclear(b[j],Rb)
    return V 

def Vr_wfn_1(molecule,bfs,P,Rb):
    """calculates electrostatic potential in point Rb
directly from wave function. Needs molecule, basis set
and density matrix. Rb is in AU."""

    a = molecule.atoms
    Rb=tuple(Rb)
    V = 0
    # nuclear contribution
    for q in range(len(molecule)):
        V+= a[q].atno / sqrt( sum( (array(a[q].pos())-Rb)**2 ) )
    # electronic contribution 
    for i in range(len(bfs)):
        for j in range(len(bfs)):
            V += P[i,j] * bfs[i].nuclear(bfs[j],Rb)
    return V 

def Vr_camms(camm,Rb):
    """calculates electrostatic potential in point Rb 
from camm distribution. camm is the object of the class 
MULT after performing the 'camms()' method on it. 
Potential is calculated in atomic units.
Remember that camm object has to contain moments in TRACELESS 
form ('makeTracelessCAMMs()' method has to be invoked before
using this function)."""

    Ra,qa,Da,Qa,Oa = camm.ReturnCAMMs()
    V=0
    for i in range(len(Ra)):
        R=Rb-Ra[i] 
        Rab=sqrt(sum(R**2,axis=0))
        V+=qa[i]/Rab
        V+=tensordot(Da[i],R,(0,0))/Rab**3 
        V+=tensordot(R,tensordot(Qa[i],R,(1,0)),(0,0))/Rab**5 
        V+=tensordot(R,tensordot(R,tensordot(Oa[i],R,(0,0)),(0,0)),(0,0))/Rab**7 
    return V
