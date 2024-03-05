"""Module containing classes for pulses. we consider pulses in the impact approximation and without the impact approximation. also pulses for tomography and interferometry."""

import os
import sys
sys.path.append("/home/martin/work/qutip/modules")
from qutip import Qobj
import numpy as np
import math

class Pulses():
    """Class for representing pulses with time delay."""

    def __init__(self, P, t, Ptype="impulse", name=None):
        """ """

        try:
            self.P     = P
            self.t     = t
            self.nP    = len(P)
            self.Ptype = Ptype
            self.name  = name
        except len(P) != len(t):
            raise Exception("Number of pulses and time delays must match. Got", len(P), "pulses and", len(t), "delays!" )


    def copy_pulse(self):
        from copy import deepcopy
        Pcopy = deepcopy(self)
        return Pcopy

    def print_pulse_info(self):
        """ """

        if self.name is not None:
            print("Pulse info for:", self.name)
        else:
            print("Pulse info for nameless pulse")
        print("Pulse type:", self.Ptype)
        if self.Ptype == "impulse":
            print("Number of pulses:", len(self.P))
            print("Pulse strengths and delays:")
            for i in range(self.nP):
                    print("Pulse", i, "P:", self.P[i], "delay:",self.t[i])


class ImpactPulses(Pulses):
    """Class for representing pulses in the impact approximation"""

    # P:    Pulse strengths
    # t:    Pulse delays
    # name: Name given to the pulses 
    def __init__(self, P, t, name=None):
        """Init method for the impact pulses"""
        super().__init__(P, t, Ptype="impact", name=Name)



class EvolutionOperator():
    """Class for representing evolution operators in the impact approximation"""

    def __init__(self, dim=2, Otype="", name=None, odd=False):
        self.dim   = dim
        self.Otype = Otype
        self.name  = name
        self.odd = odd

    def print_operator_info(self):
        print("# Evolution operator info for :", self.name)
        print("Dimensions:", self.dim)
        print("Operator type:", self.Otype)
        if self.odd:
            print("Representing odd states")
        else:
            print("Reprenting even states")

    def rename_operator(self, new_name):
        self.name = new_name

    def update_operator_type(self, new_type):
        self.Otype = new_type


class ImpulseEvolutionOperator(EvolutionOperator):
    """Class for representing pulse evolution operators in the impact approximation"""

    def __init__(self, P, dim, name=None, odd=False):
        print("Impulse operator")
        super().__init__(dim, Otype="pulse", name=name, odd=False)
        self.P  = P
        self.Up = UI(P, dim, odd)


    def update_pulse(self, P):
        self. P = P


    def update_pulse_operator(self, P):
        self.update_pulse(P)
        self.Up = UI(P, self.dim, self.odd)


    def print_pulse_operator_info(self, supress=False):
        if not supress:
            self.print_operator_info()
        print("Pulse strength:", self.P, '(au)')



class FreeEvolutionOperator(EvolutionOperator):
    """Class for representing free evolution operators in the impact approximation"""

    def __init__(self, B=1., t=0., dim=2, name=None, odd=False):
        super().__init__(dim, Otype="free", name=name, odd=odd)
        self.B  = B
        self.t  = t
        self.Uf = Uf(B, t, dim, odd)


    def update_B(self, B):
        self.B = B


    def update_time(self, t):
        self.t = t


    def update_free_operator(self, which, value):
        if which == "B":
            self.update_B(value)
        elif which == "t":
            self.update_time(value)
        self.Uf = Uf(self.B, self.t, self.dim, self.odd)


    def print_free_operator_info(self, supress=False):
        if not supress:
            self.print_operator_info()
        print("Rotational constant:", self.B, '(au)')
        print("Time delay:", self.t, '(au)')



class FullEvolutionOperator(ImpulseEvolutionOperator, FreeEvolutionOperator):
    """Class for representing combined free and pulse evolution operators in the impact approximation"""

    def __init__(self, P=0., B=1., t=0., dim=2, name=None, odd=False):
        """
        The full operator constructor
        """
        EvolutionOperator.__init__(self, dim, Otype="full", name=name, odd=odd)
        self.Up = UI(P, dim, odd)
        self.P = P
        self.Uf = Uf(B, t, dim, odd)
        self.B = B
        self.t = t
        self.U = self.Up * self.Uf
        self.Otype = "full"


    def update_full_operator(self, P=None, which=None, value=0.):
        if P is not None:
            self.update_pulse_operator(P)
        if which is not None:
            self.update_free_operator(which, value)
        self.U = self.Up * self.Uf

    def print_full_operator_info(self):
        self.print_operator_info()
        self.print_pulse_operator_info(supress=True)
        self.print_free_operator_info(supress=True)



class EvolutionOperators(ImpulseEvolutionOperator, FreeEvolutionOperator):
    """Class for representing an evolution operator with many pulses and delays """

    # Pulses: Instances of fhe Pulses class
    # B:      Rotational constant (au)
    # dim:    Dimensions of the representation
    # name:   Optianl name given to the operator
    # odd:    If True use only odd rotational states. If False, use only even states
    def __init__(self, Pulses, B=1., dim=2, name=None, odd=False):
        """Initialization method for the full impact evolution operator"""
        from qutip import Qobj, qeye
        EvolutionOperator.__init__(self, dim, Otype='Multiple pulses', name=name, odd=odd)
        self.Pulses = Pulses
        self.B = B
        self.set_full_operators()
        #Uhold = Qobj(qeye(dim))
        #for i in range(Pulses.nP):
        #    Up    = UI(Pulses.P[i], dim, odd)
        #    Ufree = Uf(B, Pulses.t[i], dim, odd)
        #    Uhold = Up * Ufree * Uhold
        #self.U = Uhold


    def set_full_operators(self):
        """Method for setting the full impact evolution operator"""
        from qutip import Qobj, qeye
        UHold = Qobj(qeye(self.dim))
        for k in range(self.Pulses.nP):
            Up    = UI(self.Pulses.P[k], self.dim, self.odd)
            Ufree = Uf(self.B, self.Pulses.t[k], self.dim, self.odd)
            UHold = Up * Ufree * UHold
        self.U = UHold

    def update_full_operators(self, P=0., Pind=None, time=0., tind=None):
        """Method for updating the pulse strengths and delay times of the full impact evolution operator"""
        icont = 0
        jcont = 0
        if Pind is not None:
            for i in Pind:
                self.Pulses.P[i] = P[icont]
                icont += 1
        if tind is not None:
            for j in tind:
                self.Pulses.t[j] = time[jcont]
                jcont += 1
        self.set_full_operators()
        



    def update_full_operators_B(self, B):
        """Method for updating the rotational constant of the full operator"""
        self.B = B
        self.set_full_operators()




# The free Hamiltonian
# @descr: The free Hamiltonian is diagonal with eigenvalues j(j+1)B
#         The pulse only connects even states with even and odd states 
#         with odd.

"""
Methods fro defining various Hamitonians and operations on the Hamiltonians based on QuTiP Qobj
"""
def H0(B, n, odd=False, full=False):
    """
    The free Hamiltonian operator
        
        Args:
            B: float
                The rotational constant

            n; int
                The dimension

            odd: Boolean
                True for odd states. False for even states

        Returns:

            H: Qobj
                The free Hamiltonian
    """
    Hmat = np.zeros((n,n))
    if full:
        for i in range(n):
            j = float(i)
            Hmat[i,i] = j * (j + 1.) * B
    else:
        jplus = 0.
        if odd:
            jplus = 1.

        for i in range(n):
            j = 2.*float(i) + jplus
            Hmat[i,i] = j * (j + 1.) * B
    H = Qobj(Hmat)
    return H

def H0_m(B:float, jmax:int, odd=False, full=False):
    """
    The free Hamiltonian operator
        
        Args:
            B: float
                The rotational constant

            jmax; int
                The max j

            odd: Boolean
                True for odd states. False for even states

        Returns:

            H: Qobj
                The free Hamiltonian
    """
    if full:
        Hmat = np.zeros(((jmax+1)**2, (jmax+1)**2))
        for j in range(jmax+1):
            jf = float(j)
            for m in range(-j,j+1):
                i = (j+1)**2 - j + m - 1
                Hmat[i,i] = jf * (jf + 1.) * B
    else:
        n = int((jmax+1)*(jmax+2)/2)
        Hmat = np.zeros((n, n))
        jplus = 0
        if odd:
            jplus = 1
            if jmax % 2 == 0:
                raise ValueError("Expected an odd integer, got", jmax) 
        else:
            if jmax % 2 != 0:
                raise ValueError("Expected an even integer, got", jmax) 
        jsub = 0
        for j in range(jplus, jmax+1, 2):
            jp = float(j)
            if j != jplus:
                jsub += 2*(j-1) + 1
                #print(j, jsub)
            for m in range(-j, j+1):
                i = (j+1)**2 - j + m - 1 - jsub - jplus
                print(j,m,i)
                Hmat[i,i] = jp * (jp + 1.) * B
    H = Qobj(Hmat)
    return H



def get_HIntMat(n, odd=False, full=False):
    """
    The cos^2(theta) - matrix in the (j,m=0) representation
    """
    from numpy import matmul

    if full:
        Hd = get_Hdip_Mat(n)
        Hmat = matmul(Hd, Hd)
        j = float(n-1)
        Hmat[-1, -1] = (2.*j + 1.)**(-1.) * ((j + 1.)**2./(2.*j + 3.) + j**2./(2.*j - 1.))

    else:
        Hmat = np.zeros((n,n))
        jplus = 0.
        if odd:
            jplus = 1.
            
        for i in range(n):
            j = 2.*float(i) + jplus
            Hmat[i,i] = (2.*j + 1.)**(-1.) * ((j + 1.)**2./(2.*j + 3.) + j**2./(2.*j - 1.))
            if i != n-1:
                Hmat[i,i+1] = (j + 1.) * (j + 2.) / ((2.*j + 1.) * (2.*j + 3.)) * math.sqrt((2*j + 1.) / (2.*j + 5.))
                Hmat[i+1,i] = Hmat[i,i+1]

    return Hmat

def get_Hdip_Mat(n):
    """
    The cos(theta) - matrix in the (j,m=0) representation
    """
    Hmat = np.zeros((n,n))
    for i in range(n-1):
        Hmat[i, i+1] = (float(i) + 1.) / math.sqrt( (2.*float(i) + 3.) * (2.*float(i) + 1.))
        Hmat[i+1, i] = Hmat[i, i+1]

    return Hmat



# The static part of the interaction Hamiltonian
# @descr: This part of the time dependent Hamiltonian will be multiplied with
#         a Gaussian pulse to form the time dependent interaction Hamiltonian 
def H1(Da, n, odd=False, full=False):
    """
    The polarizability anisotropy Hamiltonian for $m=0$
    """
    Hm = get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    return -0.25*Da*H

# The static part of the interaction Hamiltonian when using the intensity rather than the electric field
def H1_I0(Da:float, n:int, odd=False, full=False):
    from parameters import eps0, c
    Hm = get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    return -Da / (2.*eps0*c) * H

def H2(D:float, n:int):
    """
    The dipole Hamiltonian/eps for $m=0$
    Args:
    -----
        D: float
            The dipole moment
        n; int
            The dimension of the basis

    Returns:
    --------
        H: Qobj
            The dipole Hamiltonian excluding the electric field
    """

    Hm = get_Hdip_Mat(n)
    H = Qobj(Hm)
    
    return -D*H

def H_dip(eps:float, D:float, n:int):
    """
    The dipole Hamiltonian for $m=0$
    Args:
    -----
        eps: float
            The electric field strength
        D: float
            The dipole moment
        n; int
            The dimension of the basis

    Returns:
    --------
        H: Qobj
            The dipole Hamiltonian
    """

    Dmat = H2(D, n)

    return eps*Dmat




# The Gaussian pulse
def Gauss(t:float, t0:float, I0:float, sigma:float) -> float:
    """
    Gaussian pulse

        Args:
            t: float
                Time at which to evaluate the Gaussian
            t0: float
                Time for peak of Gaussian
            I0: float
                Maximum intensity (at t0)
            sigma: float
                With of Gaussian

        Returns:
            pulse: float
                The Gaussian at the specified time t
    """
    pulse = I0*math.exp(-(t-t0)**2./(2.*sigma**2.))
    return pulse

def double_Gauss(t, t0, tau, I01, I02, sigma1, sigma2):
    pulses = I01 * math.exp(-(t - t0)**2./(2.*sigma1**2.)) + I02 * math.exp(-(t- tau)**2./(2.*sigma2**2.))
    return pulses

# The Gaussian pulse for mesolver
def Gauss_me(t, args):
    return args['I0'] * math.exp(-(t-args['t0'])**2./(2.*args['sigma']**2.))

def double_Gauss_me(t, args):
    return args['I01'] * math.exp(-(t - args['t0'])**2./(2.*args['sigma']**2.)) + args['I02'] * math.exp(-(t - args['tau'])**2./(2.*args['sigma']**2.))

# Get the intensity from fluence and sigma
def getI0(P:float, sigma:float, deltaalpha:float) -> float:
    """
    Obtains the maximum pulse intensity from the pulse strength, width and polarizability anisotropy

    Args:
        P: float
            The pulse strength
        sigma: float
            the pulse width
        deltaalpha: float
            The polarizability anisotropy

    Returns:
        I0: float
            The maximum pulse intensity
    """
    from parameters import eps0, c
    I0 = 2.*eps0*c * P / (math.sqrt(2.*math.pi) * deltaalpha * sigma)
    return I0

# Get P from I0 and sigma
def getP(deltaalpha,I0,sigma):
    """
    Obtains the pulse strength from the polarizability anisotropy, maximum pulse intensity and width

    Args:
        deltaalpha: float
            The polarizability anisotropy
        I0: float
            The maximum pulse intensity
            The pulse strength
        sigma: float
            the pulse width

    Returns:
        P: float
            The pulse strength
    """

    from parameters import eps0, c
    #P = math.sqrt(2.*math.pi) * deltaalpha/(2.*eps0*c) * I0 * sigma
    P = math.sqrt(2.*math.pi) * deltaalpha/(2.*eps0*c) * I0 * sigma
    return P

# Get P from integration
def getP_int(Da, I0, t0, sigma):
    from scipy.integrate import quad
    from parameters import eps0, c
    intI, tol = quad(Gauss, -np.inf,np.inf,args=(t0, I0, sigma))
    # Where does the factor sqrt(0.5) come from??
    #return Da / (2.*eps0*c) * math.sqrt(0.5) * intI
    return Da / (2.*eps0*c) * intI

# Obtain \sigma_I (in atomic units) from T_FWHM (fs)
def sigmaFromFWHM(FWHM):
    from parameters import fs2au
    FWHM_au = FWHM * fs2au
    sigma_I = FWHM_au / math.sqrt(8.*math.log(2.))
    return sigma_I

# Impulse approximation

# The impulse Hamiltonian
def HI(P,n, odd=False, full=False):
    Hm = -P * get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    return H

# The impulse Hamiltonian mult. by -i
def HImi(P,n, odd=False, full=False):
    Hm = 1.j*P * get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    return H

# The impulse evolution operator
def UI(P,n, odd=False, full=False):
    Hm = 1.j*P * get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    U = H.expm()
    return U

# The free evolution operator
def Uf(B, t, n, odd=False, full=False):
    H = -1.j*t*H0(B,n, odd, full)
    U = H.expm()
    return U


# Composition of two evolution operators
def U2U1(U2, U1):

    try:
        return U2*U1
    except ValueError as e:
        raise Exception(e)
    except TypeError as e:
        raise Exception(e)

# Forward evolution 
# input ket
def UFW(U,x):

    try:
        if x.type == 'oper':
            return UFWO(U, x)
            #return U*x*U.dag()
        elif x.type == 'ket':
            return UFWO(U, x*x.dag())
        else:
            raise RuntimeError("Only operators or kets can be propagated in this way")
    except TypeError as e:
        raise Exception(e)

# Forward evolution
# input operator
def UFWO(U,O):
    
    return U*O*U.dag()

# Backward evolution
# input ket
def UBW(U,x):

    try:
        if x.type == 'oper':
            return UBWO(U, x)
            #return U*x*U.dag()
        elif x.type == 'ket':
            return UBWO(U, x*x.dag())
        else:
            raise RuntimeError("Only operators or kets can be propagated in this way")
    except TypeError as e:
        raise Exception(e)


# Backward evolution
# input operator
def UBWO(U,O):
    
    return U.dag()*O*U

def Proj(O1, O2):
    O = O1*O2
    
    return O.tr()


