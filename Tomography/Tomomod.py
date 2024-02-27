"""
Module containing useful methods used by tomography
"""

import math

# Matrix parameters
c2a     = 1./3.  #
c2b     = 2./(3.*math.sqrt(5.))  #
c2d     = 11./21.
c2tr    = c2a + c2d
c2diff  = c2a - c2d
disc    = c2diff**2. + 4.*c2b**2.
sqdisc  = math.sqrt(disc)
nuplus  = c2tr + sqdisc
numinus = c2tr - sqdisc
aplus   = c2a - 0.5*nuplus
nplus   = math.sqrt(aplus**2. + c2b**2.)


def Domega(P:float) -> float:
    """
    Args:
        P: float
            The pulse strength

    Returns:
        Do: float
            The angle \Delta \Omega in Eq. () in ref...
    """
    Do = -P * sqdisc
    return Do

# Star quantities

# Delta Omega
def DomegaStar():
    DoS = math.acos(-(aplus**2. - c2b**2.)**2. / (2.*aplus*c2b)**2.)
    return DoS

# Pulse strength 
def Pstar():
    DS = DomegaStar()
    PS = DS / sqdisc
    return PS

# cos(phistar)
def cosPhiS():
    Dos = DomegaStar()
    cphis = 2.*aplus*c2b * (aplus**2. - c2b**2.) * (1. - math.cos(Dos)) / nplus**4.
    return cphis


# sin(phistar)
def sinPhiS():
    Dos = DomegaStar()
    sphis = 2.*aplus*c2b * (aplus**2. + c2b**2.) * math.sin(Dos) / nplus**4.
    return sphis

# tan(phistar)
def tanPhiS():
    Dos = DomegaStar()
    tphis = (aplus**2. + c2b**2.) * math.sin(Dos) / ((aplus**2. - c2b**2.) * (1. - math.cos(Dos)))
    return tphis

def TS1(B:float, n:int) -> float:
    ps = -math.acos(cosPhiS())
    T1 = float(n)*math.pi/(3.*B) - ps/(6.*B)

    return T1

def TS2(B:float, n:int) -> float:
    T2 = math.pi/12. + TS1(B, n)
    
    return T2


def Pauli(i):
    """
    The Pauli matrices
    """
    from qutip import qeye, sigmax, sigmay, sigmaz
    try:
        if i == 0:
            Po = qeye(2)
            return Po
        elif i == 1:
            Po = sigmax()
            return Po
        elif i == 2:
            Po = sigmay()
            return Po
        elif i == 3:
            P0 = sigmaz()
            return Po
    except i != 0 or i != 1 or i !=2 or i != 3:
        raise Exception('Index', i, 'out of range for Pauli, i = 0,1,2,3')



# The standard rep
def GM(i):
    """
    The Gell-Mann matrices
    """
    from qutip import qeye, Qobj
    try:
        if i == 0:
            Go = math.sqrt(2./3.) * qeye(3)
            return Go
        elif i == 1:
            Go = Qobj(basis(3,0)*basis(3,1).dag() + basis(3,1)*basis(3,0).dag())
            return Go
        elif i == 2:
            Go = -1.j*Qobj(basis(3,0)*basis(3,1).dag() - basis(3,1)*basis(3,0).dag())
            return Go
        elif i == 3:
            Go = Qobj(basis(3,0)*basis(3,0).dag() - basis(3,1)*basis(3,1).dag())
            return Go
        elif i == 4:
            Go = Qobj(basis(3,0)*basis(3,2).dag() + basis(3,2)*basis(3,0).dag())
            return Go
        elif i == 5:
            Go = -1.j*Qobj(basis(3,0)*basis(3,2).dag() - basis(3,2)*basis(3,0).dag())
            return Go
        elif i == 6:
            Go = Qobj(basis(3,1)*basis(3,2).dag() + basis(3,2)*basis(3,1).dag())
            return Go
        elif i == 7:
            Go = -1.j*Qobj(basis(3,1)*basis(3,2).dag() - basis(3,2)*basis(3,1).dag())
            return Go
        if i == 8:
            Go = math.sqrt(1./3.) * Qobj(basis(3,0)*basis(3,0).dag() + basis(3,1)*basis(3,1).dag()) - \
                    2.*math.sqrt(1./3.) * Qobj(basis(3,2)*basis(3,2).dag())
            return Go
    except i not in [0, 1, 2, 3, 4 ,5, 6, 7, 8]:
        raise Exception('Index', i, 'out of range for Gell-Mann.')



# My rep
def GM_Martin(i):
    """
    The Gell-Mann matrices with modifications to lambda_0 -> sqrt(2/3)*lambda_0 + sqrt(1/3)lambda_8 and
    lambda_8 -> sqrt(1/3)lambda_0 - sqrt(2/3)lambda_8
    """
    from qutip import qeye, Qobj, basis
    try:
        if i == 0:
            Go = Qobj(basis(3,0)*basis(3,0).dag() + basis(3,1)*basis(3,1).dag())
            return Go
        elif i == 8:
            Go = math.sqrt(2.)*Qobj(basis(3,2)*basis(3,2).dag())
            return Go
        else:
            Go = GM(i)
            return Go

    except i != 0 or i != 1 or i != 2 or i != 3 or i != 3 or i != 4 or i != 5 or i != 6 or i !=6 or i != 7 or i != 8:
        raise Exception('Index', i, 'out of range for Gell-Mann.')


# Pauli in N-dim
def PauliN(i,N):
    """
    The 'Pauli matrices' in N dimensions. The upper left 2 x 2 corner is the Pauli matrix as specified by the index i and the rest i spadded with zeros

        Args:
            i: int
                Which Pauli matrix to obtain

            N: int
                Dimenstion of the matrix
    """
    from qutip import Qobj, basis
    try:
        if i == 0:
            Po = Qobj(basis(N,0) * basis(N,0).dag() + basis(N,1) * basis(N,1).dag())
            return Po
        elif i == 1:
            Po = Qobj(basis(N,0) * basis(N,1).dag() + basis(N,1) * basis(N,0).dag())
            return Po
        elif i == 2:
            Po = -1.j*(Qobj(basis(N,0) * basis(N,1).dag() - basis(N,1) * basis(N,0).dag()))
            return Po
        elif i == 3:
            Po = Qobj(basis(N,0) * basis(N,0).dag() - basis(N,1) * basis(N,1).dag())
            return Po
        else:
            raise Exception('Index', i, 'out of range for Pauli.')
    except N < 2:
        raise Exception('Dimenstion', i, 'too small. Expected >= 2.')

def GMN(i, N):
    """
    The 'Gell-Mann matrices' in N dimensions. The upper left 2 x 2 corner is the Pauli matrix as specified by the index i and the rest i spadded with zeros

        Args:
            i: int
                Which Gell-Mann matrix to obtain

            N: int
                Dimenstion of the matrix
    """
    from qutip import Qobj, basis
    try:
        if i == 0:
            Go = Qobj(math.sqrt(2./3.) * (basis(N,0) * basis(N,0).dag() + basis(N,1) * basis(N,1).dag()))
            return Go
        elif (i == 1) or (i == 2) or (i == 3):
            Go = PauliN(i,N)
            return Go
        elif i == 4:
            Go = Qobj(basis(N,0)*basis(N,2).dag() + basis(N,2)*basis(N,0).dag())
            return Go
        elif i == 5:
            Go = -1.j*Qobj(basis(N,0)*basis(N,2).dag() - basis(N,2)*basis(N,0).dag())
            return Go
        elif i == 6:
            Go = Qobj(basis(N,1)*basis(N,2).dag() + basis(N,2)*basis(N,1).dag())
            return Go
        elif i == 7:
            Go = -1.j*Qobj(basis(N,1)*basis(N,2).dag() - basis(N,2)*basis(N,1).dag())
            return Go
        elif i == 8:
            Go = Qobj(math.sqrt(1./3.)*(basis(N,0)*basis(N,0).dag() + basis(N,1)*basis(N,1).dag() - \
                    2.*basis(N,2)*basis(N,2).dag()))
            return Go
        else:
            raise Exception('Index', i, 'out of range for Gell-Mann.')
    except N < 3:
        raise Exception('Dimenstion', i, 'too small. Expected >= 3.')

def GM_MartinN(i, N):
    """
    The 'Gell-Mann matrices' in N dimensions. The upper left 2 x 2 corner is the Pauli matrix as specified by the index i and the rest i spadded with zeros

        Args:
            i: int
                Which Gell-Mann matrix to obtain

            N: int
                Dimenstion of the matrix
    """
    from qutip import Qobj, basis
    try:
        if i == 0:
            Go = math.sqrt(2./3.) * GMN(0, N) + math.sqrt(1./3.) * GMN(8, N)
            return Go
        elif i == 8:
            Go = math.sqrt(1./3.) * GMN(0, N) - math.sqrt(2./3.) * GMN(8, N)
            return Go
        elif (i == 1) or (i == 2) or (i == 3) or (i == 4) or (i == 5) or (i == 6) or (i == 7):
            Go = GMN(i,N)
            return Go
        else:
            raise Exception('Index', i, 'out of range for Gell-Mann.')
    except N < 3:
        raise Exception('Dimenstion', i, 'too small. Expected >= 3.')


