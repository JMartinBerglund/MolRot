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

    def __init__(self, P=0., B=1., t=0., eps=0., D=0., dim=2, a=0., name=None, odd=False, mrep=False, pol=None, alpha=0., beta=0.):
        """
        The full operator constructor
        """
        if mrep is True:
            EvolutionOperator.__init__(self, dim, Otype="fullEfield_m", name=name, odd=odd, mrep=mrep)
            if pol is None:
                self.Up    = UI(P, dim, odd, full=True, mrep=mrep)
            else:
                self.Up = UI_EFNP(P, dim, pol, alpha, beta, eps, D)
        else:
            EvolutionOperator.__init__(self, dim, Otype="fullEfield", name=name, odd=odd, mrep=mrep)
            self.Up    = UI(P, dim, odd, full=True, mrep=mrep)

        self.P     = P
        self.Uf    = UfE(B, t, eps, D, dim, a, mrep)
        self.B     = B
        self.a     = a
        self.t     = t
        self.eps   = eps
        self.alpha = alpha
        self.beta  = beta
        self.D     = D
        self.U     = U2U1(self.Up, self.Uf) #self.Up * self.UfE
        self.pol   = pol

    def update_pulse_operator(P=None, pol=None, alpha=None, beta=None):
        if P is not None:
            self.update_pulse(P)
        if alpha is not None:
            self.alpha = alpha
        if beta is not None:
            self.beta = beta
        if pol is None:
            self.Up    = UI(P, self.dim, self.odd, full=True, mrep=self.mrep)
        else:
            self.Up = UI_EFNP(P, self.dim, pol, alpha, beta, self.eps, self.D)



    def update_full_operator(self, P=None, which=None, value=0., pol=None, alpha=None, beta=None):
        if (P is not None) or (alpha is not None) or (beta is not none):
            self.update_pulse_operator(P, pol, alpha, beta)
        if which is not None:
            self.update_free_operator(which, value)
        self.U = U2U1(self.Up, self.Uf)

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


class Operator():
    """
    """

    def __init__(self, dim, mrep=False, mop=True, jstate1=0, mstate1=0, jstate2=0, mstate2=0, hermitean=False):
        """
        Construcor
        """
        from qutip import basis, qzero, ket2dm, tensor

        self.dim  = dim
        self.mrep = mrep
        self.mop  = mop

        if mrep:
            N = (dim+1)**2
            self.Op = qzero(N)
            if mop:
                j = jstate1
                if index(j,j) <= N - 1:
                    for m in range(-j, j+1):
                        self.Op += ket2dm((basis(N, index(j,m))))
                else:
                    print(index(j,j), N-1)
                    raise ValueError()
            else:
                jm = max(jstate1, jstate2)
                if index(jm, jm) <= N-1:
                    self.Op += basis(N, index(jstate1, mstate1)) * basis(N, index(jstate2, mstate2)).dag() #tensor((basis(N, index(jstate1,mstate1))), basis(N, index(jstate2,mstate2)))
                    if hermitean:
                        self.Op += self.Op.dag() 
                else:
                    raise ValueError()
                
        else:
            N = dim
            print(N)
            if mop:
                j = jstate1
                if j < N:
                    self.Op = ket2dm((basis(N, j)))
                else:
                    raise ValueError()
            else:
                jm = max(jstate1, jstate2)
                if jm < N:
                    self.Op = basis(N, jstate1) *  basis(N, jstate2).dag()
                    if hermitean:
                        self.Op += self.Op.dag() 
                else:
                    raise ValueError()



class QOpt():
    """
    Class for handling optimized quantum object preparation.
    """

    def __init__(self, qobj, qobj0, Pulses, B, a=0.):
        """
        The constructor
        """
        if qobj.shape[1] == 1:
            if (qobj.shape[0] == qobj0.shape[0]) and (qobj.shape[1] == 1):
                self.dim = qobj.shape[0]
                self.goal = qobj
                self.init = qobj0
                self.opt  = "state2state"
            else:
                raise Exception
        else:
            print("Not yet implemented")

        self.Pulses  = Pulses
        self.B       = B
        self.a       = a

        self.setEvOp()

    def setEvOp(self):
        """
        Set the initial evolution operator
        """
        from qutip import qeye
        Ui = qeye(self.dim)

        for i in range(len(self.Pulses.P)):
            Ufree = Uf(self.B, self.Pulses.t[i], int(math.sqrt(self.dim)-1), self.a, mrep=True)
            #Extend to non-parallel pulses
            UInt  = UI(self.Pulses.P[i], int(math.sqrt(self.dim)-1), odd=False, full=True, mrep=True)
            Ui = UInt * Ufree * Ui

        self.U = Ui



class StateOpt(QOpt):
    """
    Class for handling optimized state preparation.
    """

    def __init__(self, state, istate, Pulses, B, a=0.):
        """
        Constructor

            Args:
                state: Qobj (ket)
                    The target state
                istate: Qobj (ket)
                    The initial state
                Pulses: Pulses object
                    The initial guesses for the pulses
                B: float
                    The rotational constant
                a: float
                    The centrifugal distortion constant
        """
        if state.type == 'ket':
            super().__init__(state, istate, Pulses, B, a)
        else:
            raise Exception

        def update_EvOp(self, P, t, B, a, dim):
            """
            """
            from qutip import qeye
            Ui = qeye(dim)
            for i in range(len(P)):
                Ufree = Uf(B, t[i], dim, a, mrep=True)
                #Extend to non-parallel pulses
                UInt  = UI(P[i], dim, odd=False, full=True, mrep=True)
                Ui = UInt * Ufree * Ui
            U = Ui

            return U



        def opt_pulses(x, B, a, dim, state, istate):
            """
            The function to optimize
            """
            lP = len(x)/2
            P = np.zeros(lP)
            t = np.zeros(lP)
            for i in range(lP):
                P[i] = x[i]
                t[i] = x[i+lP+1]

            U = self.update_EvOp(P, t, B, a, dim)
            fstate = U * istate
            ov = fstate.overlap(state)
            return ov * ov.dag()

        def optimize(self):
            """
            Run the optimization
            """
            from scipy.optimize import minimize as mini
            lP = len(Pulses.P)
            x0 = np.zeros(2*lP)
            for i in range(lP):
                x0[i] = Pulses.P[i]
                x0[i+lP+1] = Pulses.t[i]


            res = mini(opt_pulses, x0=Pulses, args=(self.B, self.a, self.dim, self.state, self.istate))

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
                Hmat[i,i+1] = (j + 1.) * (j + 2.) / ((2.*j + 1.) * (2.*j + 3.)) * math.sqrt((2.*j + 1.) / (2.*j + 5.))
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

    
# Getting the eigenenegies and initial state based on diagonalization of the field Hamiltonian
def get_DipEigen(B, D, eps0, n):
    """
    Get the energies and the inti
    """
        
    H0_dip = H0(B, n, full=True) + H_dip(eps0, D, n) # The field Hamiltonian
            
    [Eigva, Eigve] = H0_dip.eigenstates() # Diagonalize the Hamiltonian
    return Eigva, Eigve

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


def Lorentz(x:float, x0:float, eps:float, gamma:float, offset=0.) -> float:
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

def UI_cossincos(P, n, odd=False, full=False):
    
    Hm = 1.j*P * get_cossincosm_Mat(n, odd, full)
    H = Qobj(Hm)
    U = H.expm()
    return U

def UI_cossinsin(P, n, odd=False, full=False):
    
    Hm = 1.j*P * get_cossinsinm_Mat(n, odd, full)
    H = Qobj(Hm)
    U = H.expm()
    return U


def UI_sin2cos2(P, n, odd=False, full=False):
    
    Hm = 1.j*P * 0.5 * (get_sin2diagm_Mat(n) + get_sin2cos2m_Mat(n))
    H = Qobj(Hm)
    U = H.expm()
    return U


def UI_sin2sin2(P, n, odd=False, full=False):
    
    Hm = 1.j*P * 0.5 * (get_sin2diagm_Mat(n) - get_sin2cos2m_Mat(n))
    H = Qobj(Hm)
    U = H.expm()
    return U


def UI_sin2sin2p(P, n, odd=False, full=False):
    
    Hm = 1.j*P * get_sin2sin2m_Mat(n)
    H = Qobj(Hm)
    U = H.expm()
    return U


def UI_EFNP(P:float, jmax:int, pol='z', alpha=0., beta=0., eps=0., D=0.):
    """
    The pulse propagation operator for minterferometry measurment of the local e-field with non-paralell pulses.
    """
    ca = math.cos(alpha)
    cb = math.cos(beta)
    if pol == 'z':
        Hm = 1.j*P * (ca**2. * get_Hpolm_Mat(jmax, odd=False, full=True) - \
                2.*ca*sa * get_cossincosm_Mat(jmax, odd=False, full=True) + \
                sa**2. * 0.5*(get_sin2diagm_Mat(jmax) + get_sin2cos2m_Mat(jmax, odd=False, full=True)))
    elif pol == 'x':
        sa = math.sin(alpha)
        sb = math.sin(beta)
        Hm = 1.j*P * (sa**2.*cb**2. * get_Hpolm_Mat(jmax, odd=False, full=True) + \
                2.*ca*sa*cb**2. * get_cossincosm_Mat(jmax, odd=False, full=True) - \
                2.*sa*cb*sb * get_cossinsinm_Mat(jmax, odd=False, full=True) + \
                ca**2.*cb**2. * (get_sin2_diagm_Mat(jmax) + get_sin2cos2m_Mat(jmax, odd=False, full=True)) -\
                2.*ca*cb*sb*0.5*get_sin2sin2m_Mat(jmax, odd=False, full=True) + \
                sb**2. * 0.5*(get_sin2_diagm_Mat(jmax) - get_sin2cos2m_Mat(jmax, odd=False, full=True)))

    elif pol == 'y':
        sa = math.sin(alpha)
        sb = math.sin(beta)
        Hm = 1.j*P * (sa**2.*sb**2. * get_Hpolm_Mat(jmax, odd=False, full=True) + \
                2.*ca*sa*sb**2. * get_cossincosm_Mat(jmax, odd=False, full=True) + \
                2.*sa*cb*sb * get_cossinsinm_Mat(jmax, odd=False, full=True) + \
                ca**2.*sb**2. * (get_sin2_diagm_Mat(jmax) + get_sin2cos2m_Mat(jmax, odd=False, full=True)) +\
                2.*ca*cb*sb*0.5*get_sin2sin2m_Mat(jmax, odd=False, full=True) + \
                cb**2. * 0.5*(get_sin2_diagm_Mat(jmax) - get_sin2cos2m_Mat(jmax, odd=False, full=True)))

        
    He = 1.j*H_dipm(eps, D, jmax)
    H = Qobj(Hm+He)
    U = H.expm()
    return U
                


# The free evolution operator
def Uf(B, t, n, a=0., odd=False, full=False, mrep=False):
    if mrep:
        print(mrep)
        H = -1.j*t*H0_m(B, n, a, odd, full=True)
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
            me = Gaunt_cos(j, j+1, m)    #jp1(j, m) #math.sqrt((jf + mf + 1.) * (jf - mf + 1.)/((2.*jf + 3.) * (2.*jf + 1.)))
            Hmat[i1, i2] = me
            Hmat[i2, i1] = me

    return Hmat

def get_Hsinexpplusm_Mat(jmax:int):
    """
    """
    Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
    Hmat[0,1] = Gaunt_sin(0, 0, 1, -1) #jm1mp1(1,-1)
    Hmat[3,0] = Gaunt_sin(1,1,0,0) #jp1mp1(0,0)
    # The maximum j in the representaiton separately
    for i in range(2*jmax-1):
        Hmat[index(jmax, jmax-i), index(jmax-1, jmax-1-i)] = Gaunt_sin(jmax, jmax-i, jmax-1, jmax-1-i)    #p1mp1(jmax-1, jmax-1-i)
    for j in range(1, jmax):
        jf = float(j)
        for m in range(-j, j+1):
            mf = float(m)
            i1 = index(j,m)
            i2 = index(j+1,m-1)
            Hmat[i1, i2] = Gaunt_sin(j, m, j+1, m-1) #jm1mp1(j+1, m-1)
            if (abs(m-1) <= j-1) and (j>1):
                i3 = index(j-1,m-1)
                Hmat[i1, i3] = Gaunt_sin(j,m,j-1,m-1) #jp1mp1(j-1, m-1)
    return Hmat


def get_Hsinexpminusm_Mat(jmax:int):
    """
    """

    Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
    Hmat[0,3] = Gaunt_sin(0,0,1,1) #jm1mm1(1,1)
    Hmat[1,0] = Gaunt_sin(1,-1, 0, 0) #jp1mm1(0,0)
    for i in range(2*jmax-1):
        Hmat[index(jmax, -(jmax-i)), index(jmax-1, -(jmax-1)+i)] = Gaunt_sin(jmax, -(jmax-i), jmax-1, -(jmax-1)+i)  #jp1mm1(jmax-1, -(jmax-1)+i)
    for j in range(1, jmax):
        jf = float(j)
        for m in range(-j, j+1):
            mf = float(m)
            i1 = index(j,m)
            i2 = index(j+1,m+1)
            Hmat[i1, i2] = Gaunt_sin(j+1, m+1, j, m) #jm1mm1(j+1, m+1)
            if (abs(m+1) <= j-1) and (j>1):
                i3 = index(j-1,m+1)
                Hmat[i1, i3] = Gaunt_sin(j,m,j-1,m+1) #-jp1mp1(j-1, m-1)


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

def get_cossinexpplusm_Mat(jmax:int):
    """
    """
    if jmax >= 2:
        Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
        # Take care of the j=0,1 states separately
        Hmat[0, 5] = Gaunt_cossin(0,0,2,-1) #jm1(1,0) * jm1mm1(2,1)

        for j in range(1, jmax+1):
            jf = float(j)
            for m in range(-j, j+1):
                mf = float(m)
                i1 = index(j,m)
                #f m != j:
                i2 = index(j,m-1)
                if (m > -j) or (j < jmax):
                    Hmat[i1, i2] = Gaunt_cossin(j,m,j,m-1) #jm1mm1(j+1, m+1) * jp1(j,m+1)
                if j+2 <= jmax:
                    i3 = index(j+2, m-1)
                    Hmat[i1, i3] = Gaunt_cossin(j,m,j+2,m-1) #jm1(j+1, m) * jm1mm1(j+2, m+1)
                #if (j >= 2) and (m != -j):
                if (j >= 2) and (abs(m-1) <= j-2):
                    i4 = index(j-2, m-1)
                    Hmat[i1, i4] = Gaunt_cossin(j,m,j-2,m-1) #jp1(j-1, m) * jp1mm1(j-2, m+1)

    
        return Hmat
    else:
        raise ValueError("Need jmax at least two")



def get_cossinexpminusm_Mat(jmax:int):
    """
    """
    if jmax >= 2:
        Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
        # Take care of the j=0,1 states separately
        Hmat[0, 7] = Gaunt_cossin(0,0,2,1) #jm1(1,0) * jm1mm1(2,1)

        for j in range(1, jmax+1):
            jf = float(j)
            for m in range(-j, j+1):
                mf = float(m)
                i1 = index(j,m)
                #f m != j:
                i2 = index(j,m+1)
                if (m < j) or (j < jmax):
                    Hmat[i1, i2] = Gaunt_cossin(j,m,j,m+1) #jm1mm1(j+1, m+1) * jp1(j,m+1)
                if j+2 <= jmax:
                    i3 = index(j+2, m+1)
                    Hmat[i1, i3] = Gaunt_cossin(j,m,j+2,m+1) #jm1(j+1, m) * jm1mm1(j+2, m+1)
                #if (j >= 2) and (m != -j):
                if (j >= 2) and (abs(m+1) <= j-2):
                    i4 = index(j-2, m+1)
                    Hmat[i1, i4] = Gaunt_cossin(j,m,j-2,m+1) #jp1(j-1, m) * jp1mm1(j-2, m+1)

    
        return Hmat
    else:
        raise ValueError("Need jmax at least two")



def get_cossincosm_Mat(jmax:int):
    """
    """
    Hmat = 0.5 * (get_cossinexpplusm_Mat(jmax) + get_cossinexpminusm_Mat(jmax))

    return Hmat


def get_cossinsinm_Mat(jmax:int):
    """
    """
    Hmat = -0.5j * (get_cossinexpplusm_Mat(jmax) - get_cossinexpminusm_Mat(jmax))

    return Hmat


def get_sin2m_Mat(jmax:int):
    """
    """
    if jmax >= 2:
        Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
        # Take care of the j=0,1 states separately
        Hmat[0, 0] = Gaunt_sin2(0,0,0,0) 
        Hmat[0, 4] = Gaunt_sin2(0,0,2,-2) 
        Hmat[0, 8] = Gaunt_sin2(0,0,2,2) 
        Hmat[1, 1] = Gaunt_sin2(1,-1,1,-1) 
        Hmat[2, 2] = Gaunt_sin2(1,0,1,0) 
        Hmat[3, 3] = Gaunt_sin2(1,1,1,1) 
        Hmat[1, 3] = Gaunt_sin2(1,-1, 1,1)
        Hmat[3, 1]  = Hmat[1,3]
        if jmax > 2:
            Hmat[1,9] = Gaunt_sin2(1,-1, 3,-3)
            Hmat[1,13] = Gaunt_sin2(1,-1, 3,1)
            Hmat[2,10] = Gaunt_sin2(1,0, 3,-2)
            Hmat[2,14] = Gaunt_sin2(1,0, 3,2)
            Hmat[3,11] = Gaunt_sin2(1,1, 3,-1)
            Hmat[3,15] = Gaunt_sin2(1,1, 3,3)

        for j in range(2, jmax+1):
            for m in range(-j, j+1):
                i1 = index(j,m)
                Hmat[i1,i1] = Gaunt_sin2(j,m,j,m)
                if m < j - 1:
                    i2 = index(j, m+2)
                    Hmat[i1,i2] = Gaunt_sin2(j,m,j,m+2)
                if m > -j + 1:
                    i3 = index(j, m-2)
                    Hmat[i1,i3] = Gaunt_sin2(j,m,j,m-2)
                if j < jmax - 1:
                    i4 = index(j+2,m+2)
                    Hmat[i1, i4] = Gaunt_sin2(j,m,j+2,m+2)
                    i5 = index(j+2,m-2)
                    Hmat[i1, i5] = Gaunt_sin2(j,m,j+2,m-2)
                if m < jmax - 1:
                    i6 = index(j-2, m+2)
                    Hmat[i1, i6] = Gaunt_sin2(j,m,j-2,m+2)
                if m > -jmax + 1:
                    i7 = index(j-2, m-2)
                    Hmat[i1, i7] = Gaunt_sin2(j,m,j-2,m-2)

        return Hmat
    else:
        raise ValueError("Need jmax at least 2")    


def get_sin2diagm_Mat(jmax:int):
    """
    """
    if jmax >= 2:
        Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
        # Take care of the j=0,1 states separately
        Hmat[0, 0] = Gaunt_sin2(0,0,0,0) 
        Hmat[1, 1] = Gaunt_sin2(1,-1,1,-1) 
        Hmat[2, 2] = Gaunt_sin2(1,0,1,0) 
        Hmat[3, 3] = Gaunt_sin2(1,1,1,1) 
        if jmax > 2:

            for j in range(2, jmax+1):
                for m in range(-j, j+1):
                    i1 = index(j,m)
                    Hmat[i1,i1] = Gaunt_sin2(j,m,j,m)

        return Hmat
    else:
        raise ValueError("Need jmax at least 2")    


def get_sin2exp2plusm_Mat(jmax:int):
    """
    """
    if jmax >= 2:
        Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
        # Take care of the j=0,1 states separately
        Hmat[0, 4] = Gaunt_sin2(0,0,2,-2) 
        Hmat[3, 1] = Gaunt_sin2(1,1, 1,-1)
        if jmax > 2:
            Hmat[1,9] = Gaunt_sin2(1,-1, 3,-3)
            Hmat[2,10] = Gaunt_sin2(1,0, 3,-2)
            Hmat[3,11] = Gaunt_sin2(1,1, 3,-1)

        for j in range(2, jmax+1):
            for m in range(-j, j+1):
                i1 = index(j,m)
                if m > -j + 1:
                    i2 = index(j, m-2)
                    Hmat[i1,i2] = Gaunt_sin2(j,m,j,m-2)
                if j < jmax - 1:
                    i3 = index(j+2,m-2)
                    Hmat[i1, i3] = Gaunt_sin2(j,m,j+2,m-2)
                if m > -jmax + 1:
                    i4 = index(j-2, m-2)
                    Hmat[i1, i4] = Gaunt_sin2(j,m,j-2,m-2)

        return Hmat
    else:
        raise ValueError("Need jmax at least 2")    


def get_sin2exp2minusm_Mat(jmax:int):
    """
    """
    if jmax >= 2:
        Hmat = np.zeros(((jmax + 1)**2, (jmax + 1)**2))
        # Take care of the j=0,1 states separately
        Hmat[0, 8] = Gaunt_sin2(0,0,2,2) 
        Hmat[1, 3] = Gaunt_sin2(1,-1, 1,1)
        if jmax > 2:
            Hmat[1,13] = Gaunt_sin2(1,-1, 3,1)
            Hmat[2,14] = Gaunt_sin2(1,0, 3,2)
            Hmat[3,15] = Gaunt_sin2(1,1, 3,3)

        for j in range(2, jmax+1):
            for m in range(-j, j+1):
                i1 = index(j,m)
                if m < j - 1:
                    i2 = index(j, m+2)
                    Hmat[i1,i2] = Gaunt_sin2(j,m,j,m+2)
                if j < jmax - 1:
                    i3 = index(j+2,m+2)
                    Hmat[i1, i3] = Gaunt_sin2(j,m,j+2,m+2)
                if m < jmax - 1:
                    i4 = index(j-2, m+2)
                    Hmat[i1, i4] = Gaunt_sin2(j,m,j-2,m+2)

        return Hmat
    else:
        raise ValueError("Need jmax at least 2")    

def get_sin2cos2m_Mat(jmax:int):
    """
    """
    Mat = 0.5 * (get_sin2exp2plusm_Mat(jmax) + get_sin2exp2minusm_Mat(jmax))

    return Mat


def get_sin2sin2m_Mat(jmax:int):
    """
    """
    Mat = -0.5j * (get_sin2exp2plusm_Mat(jmax) - get_sin2exp2minusm_Mat(jmax))

    return Mat
 

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

def H1_m(Da, jmax, odd=False, full=True, theta=0.):
    """
    The polarizability anisotropy Hamiltonian for $m=0$
    """
    Hm = math.cos(theta)**(2.) * get_Hpolm_Mat(jmax, odd, full, theta)
    H = Qobj(Hm)
    return -0.25*Da*H

# The static part of the interaction Hamiltonian when using the intensity rather than the electric field
def H1_I0_m(Da:float, jmax:int, odd=False, full=True, theta=0.):
    from parameters import eps0, c
    Hm = math.cos(theta)**(2.) * get_Hpolm_Mat(jmax, odd, full)
    H = Qobj(Hm)
    return -Da / (2.*eps0*c) * H

def H3_m(D:float, jmax:int, theta=0., phi=0.):
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

    Hm = math.sin(theta) * math.cos(phi) * get_Hsincosm_Mat(jmax)
    H = Qobj(Hm)
    
    return -D*H


def H4_m(D:float, jmax:int, theta=0., phi=0.):
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

    Hm = math.sint(theta) * math.sin(phi) * get_Hsinsinm_Mat(jmax)
    H = Qobj(Hm)
    
    return -D*H

"""
Analytical things related to the impact approximation
"""

def get_MatrixParams(odd=False):
    """
    Decr.
    Args:
    -----
        odd: boolean
            Whether or not to use odd states (j=1,3). Default false, i.e. we use even states (j=0,2)
    """
    if not odd:
        a, b, d = 3.**(-1.), 2.*3.**(-1.)*math.sqrt(5.)**(-1.), 11.*21.**(-1.)
    else:
        a, b, d = 3.*5.**(-1.), 2.*5.**(-1.)*math.sqrt(3./7.), 23.*45.**(-1.)

    summ = a + d
    diff = a - d
    D    = diff**2. + 4.*b**2.
    nu_p = summ + math.sqrt(D)
    nu_m = summ - math.sqrt(D)
    a_p  = a - 0.5 * nu_p
    n_p  = math.sqrt(a_p**2. + b**2.)

    return b, nu_p, nu_m, a_p, n_p
             

def MatrixElementsPulseImpact2D(P=0., odd=False):
    """
    Obtaining the matrix elements of the Pulse evolution operator in the 2D impact approximation
    Args:
    -----
        P: float
           Pulse fluence 
        Odd: boolean
             Whether or not to use odd states (j=1,3). Default false, i.e. we use even states (j=0,2)
    """
    import cmath

    b, nu_p, nu_m, a_p, n_p = get_MatrixParams(odd)

    A = (a_p**2. *  cmath.exp(1.j*0.5*P*nu_m) + b**2. * cmath.exp(1.j*0.5*P*nu_p)) * n_p**(-2.)
    B = a_p*b * (cmath.exp(1.j*0.5*P*nu_m) - cmath.exp(1.j*0.5*P*nu_p)) * n_p**(-2.)
    D = (b**2. *  cmath.exp(1.j*0.5*P*nu_m) + a_p**2. * cmath.exp(1.j*0.5*P*nu_p)) * n_p**(-2.)

    return A, B, D


def EigenField2D(B:float, D:float, eps0:float):
    """
    Obtaining the eigenenergies in of the dipole interaction in a 2-level model
    """

    V   = D * eps0 * math.sqrt(3.)**(-1.)
    l_m = B - math.sqrt(B**2. + V**2.)
    l_p = B + math.sqrt(B**2. + V**2.)

    return l_m, l_p, V


"""
Some fitting methods
"""

def least_square_fitting(data_to_fit, reference_data):
    """
    Finds the least square error of data_to_fit compared to referemce_data

    Args:
    -----
        data_to_fit: array of floats
            data to be fitted to reference data

        reference_data: array of floats
            the reference data
    """

    # Check that both arrays are of the same lengh
    if len(data_to_fit) == len(reference_data):
        #sq_err = np.dot(data_to_fit - reference_data, data_to_fit - reference_data)
        sq_err = np.mean((data_to_fit - reference_data)**2.)
        return sq_err
    else:
        raise ValueError
