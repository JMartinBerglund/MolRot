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

    def __init__(self, dim=2, Otype="", name=None, odd=False, mrep=False):
        self.dim   = dim
        self.Otype = Otype
        self.name  = name
        self.odd   = odd
        self.mrep  = mrep

    def print_operator_info(self):
        print("# Evolution operator info for :", self.name)
        print("Dimensions:", self.dim)
        print("Operator type:", self.Otype)
        if self.odd:
            print("Representing odd states")
        else:
            print("Reprenting even states")
        if self.mrep:
            print("Using non-zero m-states")
        else:
            print("Using only-zero m-states")


    def rename_operator(self, new_name):
        self.name = new_name

    def update_operator_type(self, new_type):
        self.Otype = new_type


class ImpulseEvolutionOperator(EvolutionOperator):
    """Class for representing pulse evolution operators in the impact approximation"""

    def __init__(self, P, dim, name=None, odd=False, mrep=False):
        if mrep:
            super().__init__(dim, Otype="pulse_m", name=name, odd=odd, mrep=mrep)
        else:
            super().__init__(dim, Otype="pulse", name=name, odd=odd, mrep=mrep)
        self.P  = P
        self.Up = UI(P, dim, odd, mrep)


    def update_pulse(self, P):
        self. P = P


    def update_pulse_operator(self, P):
        self.update_pulse(P)
        self.Up = UI(P, self.dim, self.odd, self.mrep)


    def print_pulse_operator_info(self, supress=False):
        if not supress:
            self.print_operator_info()
        print("Pulse strength:", self.P, '(au)')



class FreeEvolutionOperator(EvolutionOperator):
    """Class for representing free evolution operators in the impact approximation"""

    def __init__(self, B=1., t=0., dim=2, a=0., name=None, odd=False, mrep=False):
        if mrep:
            super().__init__(dim, Otype="free_m", name=name, odd=odd)
        else:
            super().__init__(dim, Otype="free", name=name, odd=odd)
        self.B  = B
        self.t  = t
        self.a  = a
        self.Uf = Uf(B, t, dim, a, odd, mrep)


    def update_B(self, B):
        self.B = B


    def update_time(self, t):
        self.t = t

    def update_a(self, a):
        self.a = a

    def update_free_operator(self, which, value):
        if which == "B":
            self.update_B(value)
        elif which == "t":
            self.update_time(value)
        elif which == "a":
            self.update_a(value)
        self.Uf = Uf(self.B, self.t, self.dim, self.a, self.odd, self.mrep)


    def print_free_operator_info(self, supress=False):
        if not supress:
            self.print_operator_info()
        print("Rotational constant:", self.B, '(au)')
        print("Centrifugal constant::", self.a, '(au)')
        print("Time delay:", self.t, '(au)')

class FreeEfieldEvolutionOperator(EvolutionOperator):
    """Class for representing free evolution operators in the impact approximation"""

    def __init__(self, B=1., t=0., eps=0., D=0., dim=2, a=0., name=None, odd=False, mrep=False):
        if mrep:
            super().__init__(dim, Otype="freeEfield_m", name=name, odd=odd)
        else:
            super().__init__(dim, Otype="freeEfield", name=name, odd=odd)
        self.B    = B
        self.t    = t
        self.a    = a
        self.eps  = eps
        self.D    = D
        self.Uf   = UfE(B, t, eps, D, dim, a, mrep)


    def update_B(self, B):
        """
        Update the rotational constant
        """
        try:
            self.B = B
        except TypeError as e:
            raise Exception(e)

    def update_a(self, a):
        """
        Update the centrifugal distortional constant
        """
        try:
            self.a = a
        except TypeError as e:
            raise Exception(e)


    def update_time(self, t):
        """
        Update the propagation time
        """
        try:
            self.t = t
        except TypeError as e:
            raise Exception(e)

    def update_eps(self, eps):
        """
        Update the electric field strength
        """
        try:
            self.eps = eps
        except TypeError as e:
            raise Exception(e)

    def update_D(self, D):
        """
        Update the dipole moment
        """
        try:
            self.D = D
        except TypeError as e:
            raise Exception(e)


    def update_free_operator(self, which, value):
        """
        Updates the free evolution operator
        """
        if which == "B":
            self.update_B(value)
        elif which == "a":
            self.update_a(value)
        elif which == "t":
            self.update_time(value)
        elif which == 'eps':
            self.update_eps(value)
        elif which == 'D':
            self.update_D(value)
        self.Uf = UfE(self.B, self.t, self.eps, self.D, self.dim, self.a, self.mrep)


    def print_free_operator_info(self, supress=False):
        if not supress:
            self.print_operator_info()
        print("Rotational constant:", self.B, '(au)')
        print("Centrifugal distortional constant:", self.a, '(au)')
        print("Time delay:", self.t, '(au)')
        print("Efield:", self.eps, '(au)')
        print("Dipole moment:", self.D, '(au)')



class FullEvolutionOperator(ImpulseEvolutionOperator, FreeEvolutionOperator):
    """Class for representing combined free and pulse evolution operators in the impact approximation"""

    def __init__(self, P=0., B=1., t=0., dim=2, a=0., name=None, odd=False, mrep=False):
        """
        The full operator constructor
        """
        if mrep:
            EvolutionOperator.__init__(self, dim, Otype="full_m", name=name, odd=odd, mrep=mrep)
        else:
            EvolutionOperator.__init__(self, dim, Otype="full", name=name, odd=odd, mrep=mrep)
        self.Up = UI(P, dim, odd, mrep)
        self.P = P
        self.Uf = Uf(B, t, dim, a, odd, mrep)
        self.B = B
        self.a = a
        self.t = t
        self.U = self.Up * self.Uf


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


class FullEfieldEvolutionOperator(ImpulseEvolutionOperator, FreeEfieldEvolutionOperator):
    """Class for representing combined free and pulse evolution operators in the impact approximation"""

    def __init__(self, P=0., B=1., t=0., eps=0., D=0., dim=2, a=0., name=None, odd=False, mrep=False):
        """
        The full operator constructor
        """
        if mrep is True:
            EvolutionOperator.__init__(self, dim, Otype="fullEfield_m", name=name, odd=odd, mrep=mrep)
        else:
            EvolutionOperator.__init__(self, dim, Otype="fullEfield", name=name, odd=odd, mrep=mrep)
        self.Up    = UI(P, dim, odd, full=True, mrep=mrep)
        self.P     = P
        self.Uf    = UfE(B, t, eps, D, dim, a, mrep)
        self.B     = B
        self.a     = a
        self.t     = t
        self.eps   = eps
        self.D     = D
        self.U     = U2U1(self.Up, self.Uf) #self.Up * self.UfE


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
    def __init__(self, Pulses, B=1., dim=2, a=0., name=None, odd=False, mrep=False):
        """Initialization method for the full impact evolution operator"""
        from qutip import Qobj, qeye
        if mrep:
            EvolutionOperator.__init__(self, dim, Otype='Multiple pulses m', name=name, odd=odd, mrep=mrep)
        else:
            EvolutionOperator.__init__(self, dim, Otype='Multiple pulses', name=name, odd=odd, mrep=mrep)
        self.Pulses = Pulses
        self.B = B
        self.a = a
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
            Up    = UI(self.Pulses.P[k], self.dim, self.odd, self.mrep)
            Ufree = Uf(self.B, self.Pulses.t[k], self.dim, self.a, self.odd, self.mrep)
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

    def update_full_operators_a(self, a):
        """Method for updating the rotational constant of the full operator"""
        self.a = a
        self.set_full_operators()


# The free Hamiltonian
# @descr: The free Hamiltonian is diagonal with eigenvalues j(j+1)B
#         The pulse only connects even states with even and odd states 
#         with odd.

"""
Methods fro defining various Hamitonians and operations on the Hamiltonians based on QuTiP Qobj
"""
def H0(B, n, a=0., odd=False, full=False):
    """
    The free Hamiltonian operator
        
        Args:
            B: float
                The rotational constant

            n: int
                The dimension

            a: float
                The rotational distorstion constant

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
            if a > 0.:
                Hmat[i,i] = j * (j + 1.) * B + j**2. * (j + 1.)**2. * a
            else:
                Hmat[i,i] = j * (j + 1.) * B
    else:
        jplus = 0.
        if odd:
            jplus = 1.

        for i in range(n):
            j = 2.*float(i) + jplus
            if a > 0.:
                Hmat[i,i] = j * (j + 1.) * B + j**2. * (j + 1.)**2. * a
            else:
                Hmat[i,i] = j * (j + 1.) * B
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


def Lorentz(x:float, x0:float, eps:float, gamma:float, offset:float) -> float:
    """
    The Lorentzian profile

    Args:
    -----
        x: float
            The argument at which to evalueate the Lorentzian
        x0: float
            The center of the Lorentzian
        eps: float
            The maximum of the Lorentizan
        gamma: float
            The widt of the Lorentzian
        offset: float
            Shift in the vertical

    Returns:
    --------
        Lorentz: float
            The Lorentzian profile
    """

    Lorentz = eps * 0.25 * gamma**2. / ((x - x0)**2. + 0.25 * gamma**2.) + offset

    return Lorentz

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
def HI(P, n, odd=False, full=False, mrep=False):
    if mrep:
        Hm = -P * get_Hpolm_Mat(n, odd, full) 
    else:
        Hm = -P * get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    return H

# The impulse Hamiltonian mult. by -i
def HImi(P, n, odd=False, full=False, mrep=False):
    if mrep:
        Hm = 1.j*P * get_Hpolm_Mat(n, odd, full)
    else:
        Hm = 1.j*P * get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    return H

# The impulse evolution operator
# Consider using HImi for this
def UI(P, n, odd=False, full=False, mrep=False):
    if mrep:
        Hm = 1.j*P * get_Hpolm_Mat(n, odd, full)
    else:
        Hm = 1.j*P * get_HIntMat(n, odd, full)
    H = Qobj(Hm)
    U = H.expm()
    return U

# The free evolution operator
def Uf(B, t, n, a=0., odd=False, full=False, mrep=False):
    if mrep:
        print(mrep)
        H = -1.j*t*H0_m(B, n, a, odd, full)
    else:
        H = -1.j*t*H0(B, n, a, odd, full)
    U = H.expm()
    return U

# The free evolution operator with constant electric field
def UfE(B, t, eps, D, n, a=0., mrep=False):
    if mrep:
        Hrot   = H0_m(B, n, a=a, odd=False, full=True)
        Hdip   = H_dipm(eps, D, n)
        H = Hrot + Hdip
    else:
        Hrot   = H0(B, n, a=a, odd=False, full=True)
        Hdip   = eps * H2(D, n)
        H = Hrot + Hdip
    #print("Hrot:", Hrot)
    #print("Hdip:", Hdip)
    #H      = Hrot + Hdip
    O      = -1.j* H * t

    U = O.expm()
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

# Thing with non-zero m

def msign(m:int) ->float:
    """
    """
    mfact = float(m) + 0.5*float(m)*(1. - np.sign(m))
    msign = (-1.)**mfact

    return msign

def index(j:int, m:int) -> int:
    """
    Retruns the index of the matrix or vector

    Args:
    -----
        j: int
            The j quantum number
        m: int
            The m quantum number

    Returns:
    --------
        i: int
            The index
    """

    i = j * (j + 1) + m

    return i

def jp1(j:int, m:int) -> float:
    """
    """
    jf = float(j)
    mf = float(m)

    if abs(m) <= j:
        me = (jf + mf + 1.) * (jf - mf + 1.) / ((2.*jf + 3.) * (2.*jf + 1.))
        return math.sqrt(me)
    else:
        raise ValueError("Absolute value of m cannot exeed j")

def jm1(j:int, m:int) -> float:
    """
    """
    jf = float(j)
    mf = float(m)

    if abs(m) <= j:
        me = (jf + mf) * (jf - mf) / ((2.*jf + 1.) * (2.*jf - 1.))
        return math.sqrt(me)
    else:
        raise ValueError("Absolute value of m cannot exeed j")



def jp1mp1(j:int, m:int) -> float:
    """
    """
    jf = float(j)
    mf = float(m)
    mesign = msign(m) * msign(m+1)
    if abs(m) <= j:
        me = (jf + mf + 2.) * (jf + mf + 1.) /((2.*jf + 3.) * (2.*jf + 1.))
        return mesign * math.sqrt(me)
    else:
        raise ValueError("Absolute value of m cannot exeed j")

def jp1mm1(j:int, m:int) -> float:
    """
    """
    jf = float(j)
    mf = float(m)
    mesign = msign(m) * msign(m-1)
    if abs(m) <= j:
        me = (jf - mf + 2.) * (jf - mf + 1.) /((2.*jf + 3.) * (2.*jf + 1.))
        return mesign * math.sqrt(me)
    else:
        raise ValueError("Absolute value of m cannot exeed j")


def jm1mp1(j:int, m:int) -> float:
    """
    """
    jf = float(j)
    mf = float(m)
    mesign = msign(m) * msign(m+1)
    if abs(m) <= j:
        me = (jf - mf) * (jf - mf - 1.) /((2.*jf + 1.) * (2.*jf - 1.))
        return mesign * math.sqrt(me)
    else:
        raise ValueError("Absolute value of m cannot exeed j")

def jm1mm1(j:int, m:int) -> float:
    """
    """
    jf = float(j)
    mf = float(m)
    mesign = msign(m) * msign(m-1)
    if abs(m) <= j:
        me = (jf + mf) * (jf + mf - 1.) /((2.*jf + 1.) * (2.*jf - 1.))
        return mesign * math.sqrt(me)
    else:
        raise ValueError("Absolute value of m cannot exeed j")


def Ysign(m:int) -> float:
    """
    Returns the sign of the spherical harmonic Y_{j,m}

    Args:
        m: int
            The magnetic quantum number of Y

    Returns:
        msign: float
            The sign of the spherical harmonic (1, -1)
    """

    msign = (-1.)**(m + 0.5*m*(1. - math.copysign(1,m)))

    return msign

def Gaunt_cos(j1:int, j2:int, m:int) -> float:
    """
    """
    from sympy.physics.wigner import gaunt
    if abs(j1 - j2) == 1:
        pref = 2. * math.sqrt(math.pi/3.)
        if m >= 0:
            Y1sign = Ysign(m)
        else:
            Y1sign = Ysign(-m)
        me = pref * Y1sign * float(gaunt(j1, 1, j2, -m, 0, m))
        return me
    else: # Can't get the except to excecute!!
        raise ValueError("The selection rule for cos operator requires that \Delta j = pm 1")


def Gaunt_cos2(j1:int, j2:int, m:int) -> float:
    """
    """
    from sympy.physics.wigner import gaunt
    if (abs(j1 - j2) == 2) or (j1 == j2):
        pref1 = 4./3. * math.sqrt(math.pi/5.)
        pref2 = 2./3. * math.sqrt(math.pi)
        if m >= 0:
            Y1sign = Ysign(m)
        else:
            Y1sign = Ysign(-m)
        me = Y1sign * (pref1 * float(gaunt(j1, 2, j2, -m, 0, m)) + pref2 * float(gaunt(j1, 0, j2, -m, 0, m)))
        return me
    else: # Can't get the except to excecute!!
        raise ValueError("The selection rule for cos^2 operator requires that \Delta j = 0, pm 2")


def Gaunt_cos2_20(j1:int, j2:int, m:int) -> float:
    """
    """
    from sympy.physics.wigner import gaunt
    if abs(j1-j2) == 2:
        pref = 4./3. * math.sqrt(math.pi/5.)
        if m >= 0:
            Y1sign = Ysign(m)
        else:
            Y1sign = Ysign(m)
        me = Y1sign * pref * float(gaunt(j1, 2, j2, -m, 0, m))
        return me
    else: # Can't get the except to excecute!!
        raise ValueError("The selection rule dictates \Delta j = pm 2")



def Gaunt_cos2_00(j1:int, j2:int, m:int) -> float:
    """
    """
    from sympy.physics.wigner import gaunt
    if j1 == j2:
        pref = 2./3. * math.sqrt(math.pi)
        if m >= 0:
            Y1sign = Ysign(m)
        else:
            Y1sign = Ysign(-m)
        me = Y1sign * pref * float(gaunt(j1, 0, j2, -m, 0, m))
        return me
    else: # Can't get the except to excecute!!
        raise ValueError("The selection rule dictates \Delta j = 0")


def Gaunt_sin(j1:int, m1:int, j2:int, m2:int) -> float:
    """
    """
    from sympy.physics.wigner import gaunt
    if (abs(m1 - m2) == 1) and (abs(j1-j2) == 1):
        m_min = min(m1,m2)
        m_max = max(m1,m2)
        prefact = 2.*math.sqrt(2.*np.pi/3.)

        if m_max <= 0:
            if m1 > m2:
                me = Ysign(-m1) * prefact * (-1.) * float(gaunt(j1,1,j2,-m1,1,m2))
            else:
                me = Ysign(-m1) * prefact * float(gaunt(j1,1,j2,-m1,-1,m2))
        elif m_min <= 0:
            if m1 > m2:
                me = Ysign(m1) * prefact * (-1.) * float(gaunt(j1,1,j2,-m1,1,m2))
            else:
                me = Ysign(m1) * prefact * float(gaunt(j1,1,j2,-m1,-1,m2))
        else:
            if m1 > m2:
                me = Ysign(m1) * prefact * (-1.) * float(gaunt(j1,1,j2,-m1,1,m2))
            else:
                me = Ysign(m1) * prefact * float(gaunt(j1,1,j2,-m1,-1,m2))
        return me
                                                                                                                                                                                                            
    elif abs(m1-m2) != 1:
        raise ValueError('Violating m  selection rule')
    else:
        raise ValueError('Violating j selection rule')

def Gaunt_cossin(j1:int, m1:int, j2:int, m2:int) -> float:
    """
    """
    from sympy.physics.wigner import gaunt
    if (abs(m1-m2) == 1) and ((abs(j1-j2) == 2) or (j1 == j2)):
        m_min = min(m1,m2)
        m_max = max(m1,m2)
        prefact = 2.*math.sqrt(2.*np.pi/15.)

        if m_max <= 0:
            if m1 > m2:
                me = Ysign(-m1) * prefact * (-1.) * float(gaunt(j1,2,j2,-m1,1,m2))
            else:
                me = Ysign(-m1) * prefact * float(gaunt(j1,2,j2,-m1,-1,m2))
        elif m_min <= 0:
            if m1 > m2:
                me = Ysign(m1) * prefact * (-1.) * float(gaunt(j1,2,j2,-m1,1,m2))
            else:
                me = Ysign(m1) * prefact * float(gaunt(j1,2,j2,-m1,-1,m2))
        else:
            if m1 > m2:
                me = Ysign(m1) * prefact * (-1.) * float(gaunt(j1,2,j2,-m1,1,m2))
            else:
                me = Ysign(m1) * prefact * float(gaunt(j1,2,j2,-m1,-1,m2))
        return me
                                                                                                                                                                                                            
    elif abs(m1-m2) != 1:
        raise ValueError('Violating m selection rule')
    else:
        raise ValueError('Violating j selection rule')


def Gaunt_sin2(j1:int, m1:int, j2:int, m2:int) -> float:
    """
    """
    from sympy.physics.wigner import gaunt
    if (abs(m1-m2) == 2) and ((abs(j1-j2) == 2) or (j1 == j2)):
        m_min = min(m1,m2)
        m_max = max(m1,m2)
        prefact = 4.*math.sqrt(2.*np.pi/15.)

        if m_max <= 0:
            if m1 > m2:
                me = Ysign(-m1) * prefact * float(gaunt(j1,2,j2,-m1,2,m2))
            else:
                me = Ysign(-m1) * prefact * float(gaunt(j1,2,j2,-m1,-2,m2))
        elif m_min <= 0:
            if m1 > m2:
                me = Ysign(m1) * prefact  * float(gaunt(j1,2,j2,-m1,2,m2))
            else:
                me = Ysign(-m1) * prefact * float(gaunt(j1,2,j2,-m1,-2,m2))
        else:
            if m1 > m2:                                                                
                me = Ysign(m1) * prefact * float(gaunt(j1,2,j2,-m1,2,m2))
            else:
                me = Ysign(m1) * prefact * float(gaunt(j1,2,j2,-m1,-2,m2))
        return me

    #elif (abs(m1-m2) == 2) and (j1 == j2):
        

    elif (m1 == m1) and (j1 == j2):
        me = 1. - Gaunt_cos2(j1, j2, m1)
        return me


    #elif (m1 == m2) and (abs(j1-j2) == 2):
    #    prefact = -4./3.*math.sqrt(np.pi/5.)

    #    if m1 <= 0:
    #        me = Ysign(-m1) * prefact * float(gaunt(j1,2,j2,-m1,0,m2))
    #    else:
    #        me = Ysign(m1) * prefact  * float(gaunt(j1,2,j2,-m1,0,m2))

    #    return me

    #elif (m1 == m2) and (j1 == j2):
        #me = 1. - Gaunt_cos2(j1, j2, m1)
        #prefact = 4./3.*math.sqrt(np.pi)

        #if m1 <= 0:
        #    me = Ysign(-m1) * prefact * float(gaunt(j1,2,j2,-m1,2,m2))
        #else:
        #    me = Ysign(m1) * prefact  * float(gaunt(j1,2,j2,-m1,2,m2))
    
        #return me
                                                                                                                                                                                                            
    else:
        raise ValueError('Violating j or m selection rules')



def H0_m(B:float, jmax:int, a=0., odd=False, full=False):
    """
    The free Hamiltonian operator
        
        Args:
            B: float
                The rotational constant

            jmax; int
                The max j

            a: float
                The rotational disttorsion constant

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
                if a > 0.:
                    Hmat[i,i] = jf * (jf + 1.) * B + jf**2. * (jf + 1.)**2. * a
                else:
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
                if a > 0.:
                    Hmat[i,i] = jf * (jf + 1.) * B + jf**2. * (jf + 1.)**2. * a
                else:
                    Hmat[i,i] = jp * (jp + 1.) * B
    H = Qobj(Hmat)
    
    return H

def get_Hdipm_Mat(jmax:int):
    """
    """

    Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
    for j in range(jmax):
        jf = float(j)
        for m in range(-j, j+1):
            mf = float(m)
            i1 = index(j,m)
            i2 = index(j+1, m)
            me = jp1(j, m) #math.sqrt((jf + mf + 1.) * (jf - mf + 1.)/((2.*jf + 3.) * (2.*jf + 1.)))
            Hmat[i1, i2] = me
            Hmat[i2, i1] = me

    return Hmat

def get_Hsinexpplusm_Mat(jmax:int):
    """
    """

    Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
    Hmat[0,1] = jm1mp1(1,-1)
    Hmat[3,0] = jp1mp1(0,0)
    for i in range(jmax+1):
        Hmat[index(jmax, jmax-i), index(jmax-1, jmax-1-i)] = jp1mp1(jmax-1, jmax-1-i)
    for j in range(1, jmax):
        jf = float(j)
        for m in range(-j, j+1):
            mf = float(m)
            i1 = index(j,m)
            i2 = index(j+1,m-1)
            Hmat[i1, i2] = jm1mp1(j+1, m-1)
            if (abs(m-1) <= j-1) and (j>1):
                i3 = index(j-1,m-1)
                Hmat[i1, i3] = jp1mp1(j-1, m-1)
    return Hmat


def get_Hsinexpminusm_Mat(jmax:int):
    """
    """

    Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
    Hmat[0,3] = jm1mm1(1,1)
    Hmat[1,0] = jp1mm1(0,0)
    for i in range(jmax+1):
        Hmat[index(jmax, -(jmax-i)), index(jmax-1, -(jmax-1)+i)] = jp1mm1(jmax-1, -(jmax-1)+i)
    for j in range(1, jmax):
        jf = float(j)
        for m in range(-j, j+1):
            mf = float(m)
            i1 = index(j,m)
            i2 = index(j+1,m+1)
            Hmat[i1, i2] = jm1mm1(j+1, m+1)
            #if (abs(m-1) <= j-1) and (j>1):
            #    i3 = index(j-1,m-1)
            #    Hmat[i1, i3] = -jp1mp1(j-1, m-1)


    return Hmat

def get_Hsincosm_Mat(jmax:int):
    """
    """
    Hmat = 0.5 * (get_Hsinexpplusm_Mat(jmax) + get_Hsinexpminusm_Mat(jmax))

    return Hmat
    
def get_Hsinsinm_Mat(jmax:int):
    """
    """
    Hmat = -0.5j * (get_Hsinexpplusm_Mat(jmax) - get_Hsinexpminusm_Mat(jmax))

    return Hmat

def get_cossinexpminusm_Mat(jmax:int):
    """
    """
    if jmax >= 2:
        Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
        # Take care of the j=0,1 states separately
        Hmat[0, 7] = jm1(1,0) * jm1mm1(2,1)
        print(jm1mm1(2,0) * jp1(1,0))
        print(jm1mm1(2,1) * jp1(1,1))

        for j in range(1, jmax+1):
            print(j)
            jf = float(j)
            for m in range(-j, j):
                mf = float(m)
                i1 = index(j,m)
                #f m != j:
                i2 = index(j,m+1)
                Hmat[i1, i2] = jm1mm1(j+1, m+1) * jp1(j,m+1)
                print("1", j,m, i1, i2)
                if j+2 <= jmax:
                    i3 = index(j+2, m+1)
                    Hmat[i1, i3] = jm1(j+1, m) * jm1mm1(j+2, m+1)
                    print("2", j, m, i1, i3)
                #if (j >= 2) and (m != -j):
                if (j >= 2) and (abs(m+1) <= j-2):
                    i4 = index(j-2, m+1)
                    Hmat[i1, i4] = jp1(j-1, m) * jp1mm1(j-2, m+1)
                    print("3", j, m, i1, i4)

    
        return Hmat
    else:
        raise ValueError("Need jmax at least two")



def get_cossincos_m_Mat(jmax:int):
    """
    """
    from numpy import matmul

    Hmat = matmul(get_Hdipm_Mat(jmax), get_Hsincosm_Mat(jmax))
    i1 = index(jmax, jmax)
    i2 = index(jmax, jmax-1)
    i3 = index(jmax, -jmax)
    i4 = index(jmax, -jmax + 1)
    Hmat[i1, i2] = jm1(jmax, jmax-1) * jp1mp1(jmax-1, jmax-1)
    Hmat[i3, i4] = jm1(jmax, -jmax+1) * jp1mm1(jmax-1, -jmax+1)
    # WARNING NEED TO CORRECT THE MATRIX ELEMENTS BELONINGING TO JMAX
    print("WARNING, MATRIX ELEMENTS FOR JMAX ARE NOT CORRECT!!")

    return Hmat
 
def get_cossinsin_m_Mat(jmax:int):
    """
    """
    from numpy import matmul

    Hmat = matmul(get_Hdipm_Mat(jmax), get_Hsinsinm_Mat(jmax))

    # WARNING NEED TO CORRECT THE MATRIX ELEMENTS BELONINGING TO JMAX
    print("WARNING, MATRIX ELEMENTS FOR JMAX ARE NOT CORRECT!!")

    return Hmat
 
def H2_m(D:float, jmax:int):
    """
    The dipole Hamiltonian/eps for $m=0$
    Args:
    -----
        D: float
            The dipole moment
        jmax: int
            The maximum j for the representation

    Returns:
    --------
        H: Qobj
            The dipole Hamiltonian excluding the electric field
    """

    Hm = get_Hdipm_Mat(jmax)
    H = Qobj(Hm)
    
    return -D*H

def H_dipm(eps:float, D:float, jmax:int):
    """
    The dipole Hamiltonian for $m=0$
    Args:
    -----
        eps: float
            The electric field strength
        D: float
            The dipole moment
        jmax int
            The maximum j for the representation

    Returns:
    --------
        H: Qobj
            The dipole Hamiltonian
    """
    Dmat = H2_m(D, jmax)

    return eps*Dmat


def get_Hpolm_Mat(jmax:int, odd=False, full=True):
    """
    Obtain the cos^2(theta) matrix in the (j,m)-representation
    """
    from numpy import matmul
    
    # ADD POSSIBILITY FOR EVEN AND ODD REPRESENTATION, NOT JUST THE FULL ONE
    Hm = get_Hdipm_Mat(jmax)
    Hmat = matmul(Hm, Hm)
    jm = float(jmax)
    i1 = index(jmax, jmax)
    i2 = index(jmax, -jmax)
    Hmat[i1, i1] = (2.*jm + 1.)/((2.*jm + 3.) * (2.*jm + 1.))
    Hmat[i2, i2] = Hmat[i1, i1] 

    return Hmat

def H1_m(Da, jmax, odd=False, full=True):
    """
    The polarizability anisotropy Hamiltonian for $m=0$
    """
    Hm = get_Hpolm_Mat(jmax, odd, full)
    H = Qobj(Hm)
    return -0.25*Da*H

# The static part of the interaction Hamiltonian when using the intensity rather than the electric field
def H1_I0_m(Da:float, jmax:int, odd=False, full=True):
    from parameters import eps0, c
    Hm = get_Hpolm_Mat(jmax, odd, full)
    H = Qobj(Hm)
    return -Da / (2.*eps0*c) * H


