"""Module containing tomography classes and methods in the impulse approximation and without the approximation"""

#!/usr/bin/env python3
from textwrap import dedent
import os
import sys
sys.path.append("/home/martin/Work/quantumLibWorkP3")
import numpy as np
import math
import cmath
from qutip import *


# Parameters (these belong in paramters.py)
eps0 = (4.*math.pi)**(-1.) # vacuum permittivitiy (atomic units)
c    = 137.03599  # speed of light (atomic units)
fs2au = 41.341373 # converting femtoseconds to atomic time

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

class AnalyticalTomography():
    """
    A class for analytical tomography for a 2 level basis
    """

    def __init__(self, P=0., B=0., t=[], istate=None, Name=None) -> None:
        """
        The constructor

        Args:
        -----
            P: float
                The pulse strength
        
            B: float
                The rotational constant

            t: array of floats
                Time grid

            istate: Qobj
                The quantum object to be measured

            Name: str
                The name or the instance


        Parameters:
        -----------
            DO: float
                Delta Omega
         """
        from Utility import Uf, UI, U2U1 
        self.P = P
        self.B = B
        self.DO = Domega(P) #-P * sqdisc
        self.Stokes = [Stokesi(0, P), Stokesi(1, P), Stokesi(2, P), Stokesi(3, P)]
        if istate is None:
            from qutip import basis
            self.istate = basis(2,0)
        else:
            self.istate = istate

        # Star quantities
        self.Pstar = Pstar() #2.6726819190994378
        self.DOstar = Domega(Pstar()) #-self.Pstar * self.sqdisc
        self.T1Star = TS1(B, 0)
        self.T2Star = TS2(B, 0)
        self.Us = [
                U2U1(UI(self.Pstar, 2), Uf(self.B, self.T1Star, 2)),
                U2U1(UI(self.Pstar, 2), Uf(self.B, self.T2Star, 2)), 
                UI(0., 2)
        ]
        self.Stokes = np.zeros(4)
        self.mstate = None

        if Name is None:
            self.Name = ''
        else:
            self.Name = Name

    def set_P(self, P):
        """
        """
        try:
            self.P = P
            self.DO = Domega(P) #-P * self.sqdisc
            self.Stokes = [Stokesi(0, P), Stokesi(1, P), Stokesi(2, P), Stokesi(3, P)]
        
        except P is None:
            raise TypeError("P cannot be None")
        
        except P < 0.:
            raise ValueError("P must be positive, got", P)

    def print_Stokes(self):
        """
        """
        for i in range(4):
            print("Stokes parameter",i,":",self.Stokes[i])

    def print_iStokes(self):
        """
        """
        from Utility import Proj
        for i in range(4):
            print("Initial Stokes parameter",i,":",Proj(Pauli(i), self.istate * self.istate.dag()))


    def StokesP(self, P):
        """

        Args:
            P: array of floats
                The values of P for which to calculatet the Stokes paramters
        """
        lP = len(P)

        for i in range(lP):
            DO = Domega(P[i])
            print(i, P[i], Stokesi(0, P[i]), Stokesi(1, P[i]), Stokesi(2, P[i]), Stokesi(3, P[i]))

    def plot_StokesP(self, P):
        """
        """
        import matplotlib.pyplot as plt
        lP = len(P)
        S0 = np.ones(lP)
        S1 = np.zeros(lP)
        S2 = np.zeros(lP)
        S3 = np.zeros(lP)

        for i in range(lP):
            S1[i] = Stokesi(1, P[i])
            S2[i] = Stokesi(2, P[i])
            S3[i] = Stokesi(3, P[i])

        plt.figure()
        plt.xlabel("P (atomic units)")
        plt.ylabel("Projection on Stokes matrix")
        plt.plot(P, S0, label='S0')
        plt.plot(P, S1, '--', label='S1')
        plt.plot(P, S2, label='S2')
        plt.plot(P, S3, label='S3')
        plt.legend(loc='upper right')
        plt.show()

    def tomo_run(self):
        from qutip import basis
        from Utility import Proj
        Op = basis(2, 0) * basis(2,0).dag()
        self.Stokes[0] = 1.
        for i in range(3):
            fstate = self.Us[i] * self.istate
            fspop = Proj(Op, fstate * fstate.dag()).real
            self.Stokes[i+1] = 2. * fspop - 1.        



class Tomography():
    """
    The main tomography class
    """

    def __init__(self, Pulsepara=None, Molpara=None, istate=None, dim=3, Methods=None, Lagrange=None, Name=None, Ptype="impulse") -> None:
        """Initialization routine
        
        Args:
            Pulsepara: Parafile
                Containing the pulse fluencies and time delays

            B: float
                The rotational constant

            dim: int 
                The dimension of the basis used

            Methods: dict 
                Which method is used to obtain the tomographic pulses

            Lagrange: dict 
                The Lagrange multipliers used for optimization

            Name: string
                Optional name for the tomography run
        """
        from Utility import Pulses, GaussPulses
        from parameters import Molecule
        self.dim      = dim
        self.istate   = istate
        self.tstate   = None
        self.Methods  = Methods
        self.Lagrange = Lagrange
        self.Name     = Name
        self.Ptype    = Ptype
        self.pops     = None
        self.Ops      = None
        
        if Pulsepara is not None:
            #Create an instance of the Pulses class. Each P and \tau pair belongs to one pulse sequence in the lab
            if "Pulses" in Pulsepara:
                nP = len(Pulsepara['Pulses'])
                self.Pulses =[]
                if Ptype == "impulse":
                    for i in range(nP):
                        Para = Pulsepara['Pulses'][i]
                        self.Pulses.append(Pulses(Para, Ptype))
                elif Ptype == "Gauss":
                    for i in range(nP):
                        Para = Pulsepara['Pulses'][i]
                        self.Pulses.append(GaussPulses(Para, Ptype))
            else:
                if Ptype == "impulse":
                    self.Pulses = Pulses(Pulsepara, Ptype)
                elif Ptype == "Gauss":
                    self.Pulses = GaussPulses(Pulsepara, Name)
            #self.set_pulses(Pulses)
        else:
            self.Pulses = None
            self.Us = None

        if Molpara is None:
            self.Molecule = Molecule({'B': 1., 'Da': 1., 'Mol': "Custom"})
        else:
            if isinstance(Molpara, dict):
                if "B" not in Molpara:
                    Molpara.update({"B": 1.})
                if "Da" not in Molpara:
                    Molpara.update({"Da": 1.})
                if 'Mol' not in Molpara:
                    Molpara.update({'Mol': "Custom"})
            self.Molecule = Molecule(Molpara)
        self.B = self.Molecule.B # Remove this at some point, B is already in Molecule

    def print_info(self, full_mol_info=False):
        """
        """
        if self.Name is None:
            print("Printing info on nameless {}Tomography instance".format(self.Ptype))
        else:
            print("Printing info for {}Tomography instance called".format(self.Ptype, self.Name))
        print()
        print("Molecule info:")
        print("--------------")
        self.Molecule.print_info(full_mol_info)
        print()
        print("Pulse info:")
        print("-----------")
        if isinstance(self.Pulses, list):
            for i in range(len(self.Pulses)):
                self.Pulses[i].print_pulse_info()
                print()
        else:
            self.Pulses.print_pulse_info()
            print()


    @staticmethod
    def get_Stokes_projs(self, rho, dimStokes):
        """Obtain the Stokes paramters, {S_i} from the density matrix
        Args:
            rho : QuTip Qobj
                The density matrix
        """
        if dimStokes is not None:
            if dimStokes <= rho.dims**2:
                Stokes = np.zeros(dimStokes)
                for i in range(dimStokes):
                    Op = PauliN(i, self.dim)
                    Stokes[i] = Proj(rho, Op)
                return Stokes
            else:
                raise ValueError("The number of Stokes parameters cannot exceed the square of the dimension of the quantum object")
        else:
            raise RuntimeError("Need to know the number of Stokes paramters, none given")


    def get_Stokes_projections(self, rho, dimStokes):
        """Obtain the Stokes paramters, {S_i} from the density matrix
        Args:
            rho : QuTip Qobj
                The density matrix
        """
        from Utility import Proj
        if dimStokes is not None:
            if dimStokes <= rho.dims[0][0]**2:
                Stokes = np.zeros(dimStokes)
                for i in range(dimStokes):
                    Op = PauliN(i, self.dim)
                    Stokes[i] = Proj(rho, Op)
                return Stokes
            else:
                raise ValueError("The number of Stokes parameters cannot exceed the square of the dimension of the quantum object")
        else:
            raise RuntimeError("Need to know the number of Stokes paramters, none given")





class ImpactTomography(Tomography):
    """
    Implements tomography in the impact approximation
    """

    def __init__(self, Pulsepara=None, Molpara=None, istate=None, dim=3, Methods=None, Lagrange=None, Name=None) -> None:
        """Initialization routine
        
        Args:
            Pulsepara: Parafile
                Containing the pulse fluencies and time delays

            Molpara: Parafile
                The molecular parameters (e.g. rotational constant B, polarizability anisotropy Da)

            istate: Qobj:int
                The iniital state of the molecule

            dim: int 
                The dimension of the basis used

            Methods: dict 
                Which method is used to obtain the tomographic pulses

            Lagrange: dict 
                The Lagrange multipliers used for optimization

            Name: string
                Optional name for the tomography run
        """
        try:
            # Call the main initi routine
            super().__init__(Pulsepara, Molpara, istate, dim, Methods, Lagrange, Name, Ptype="impulse")
            # Next we set the evolution operators corresponding to each pulse sequence
            if Pulsepara is not None:
                self.set_EvolutionOperators()
        except KeyError as ke:
            print("KeyError: {} occurred in ImpactTomography".format(ke))
        except TypeError as te:
            print("TypeError: {} occurred in ImpactTomography".format(te))
        except ValueError as ve:
            print("ValueError: {} occurred in ImpactTomography".format(ve))
        except Exception as e:
            print("Exception: {} occurred in ImpactTomography".format(e))


    def set_pulses(self, Pulses) -> None:
        """Sets the pulses for the tomography instance
        
        Args: 
            Pulses : Pulsepara containing the pulses strengths and time delays

        """
        import Utility as Ut
        nP  = len(Pulses['Ps'])
        self.Pulses  = []
        for i in range(nP):
            Para = {'Ps': Pulses['Ps'][i], 'taus': Pulses['taus'][i]}
            #Pulse = Ut.Pulses(Pulses[i].P, Pulses[i].t)
            Pulse = Ut.Pulses(Para, "impulse")
            self.Pulses.append(Pulse)

    def set_EvolutionOperators(self) -> None:
        """Sets the three evolution operators for the tomography scheme"""
        import Utility as Ut
        #nP  = len(self.Pulses.Ps)
        self.Us  = []
        if isinstance(self.Pulses, list):
            lP = len(self.Pulses)
            for i in range(lP):
                lPs = len(self.Pulses[i].Ps)
                for j in range(lPs):
                    Para = {'Ps': [self.Pulses[i].Ps[j]], 'taus': [self.Pulses[i].taus[j]]}
                    Pulse = Ut.Pulses(Para, "impulse")
                    self.Us.append(Ut.EvolutionOperators(Pulse, self.B, dim=self.dim))
        else:
            for i in range(self.Pulses.nP):
                # Create Para-files for the EvolutionOperators class
                Para = {'Ps': [self.Pulses.Ps[i]], 'taus': [self.Pulses.taus[i]]}
                Pulse = Ut.Pulses(Para, "impulse")
                #self.Us.append(Ut.EvolutionOperators(self.Pulses.Ps[i], self.B, dim=self.dim))
                self.Us.append(Ut.EvolutionOperators(Pulse, self.B, dim=self.dim))


    @staticmethod
    def pop_run(self, i:int) -> float:
        """Running a measurment for pulse i
        
        Args:
            i : int
                Which pulse to chose
        
        Returns:
            pop0 : float
                The final ground state population
        """
        in_st  = self.istate
        U      = self.Us[i].U
        fin_st = U * in_st
        ov     = fin_st.overlap(basis(self.dim,0))
        pop0   = ov * np.conj(ov)

        return pop0.real

    def pop_runs(self) -> None:
        """Obtain the ground state populations from forward propagations with various pulses."""

        npop = len(self.Us)
        pop0s = np.zeros(npop)
        for i in range(npop):
            pop0s[i] = self.pop_run(i)

        self.pop0s = pop0s
    
    def BW_prop(self, i:int, Op=None):
        """
        """
        from Utility import UBWO
        if Op is None:
            Op = PauliN(3, self.dim)
        U = self.Us[i].U
        Op_BW = UBWO(U, Op)

        return Op_BW

    def BW_prop_all(self, Op):
        """
        """
        from Utility import UBWO, U2U1
        from qutip import qeye
        if Op is None:
            Op = PauliN(3, self.dim)
        U = qeye(self.dim)
        for i in range(len(self.Us)):
            U = U2U1(self.Us[i].U, U)
            #U = U2U1(U, self.Us[i].U)
        Op_BW = UBWO(U, Op)

        return Op_BW



    def print_pulses(self, ind) -> None:
        """
        """
        import matplotlib.pyplot as plt

        try:
            plt.figure()
            if isinstance(self.Pulses, list):
                nPu = len(self.Pulses[ind].taus)
                if nPu == 1:
                    print("One pulse list")
                    tau = self.Pulses[ind].taus[0]
                    P   = self.Pulses[ind].Ps[0] 
                    xmin = tau - 1.  
                    xmax = tau + 1.  
                    plt.xlim(xmin, xmax)
                    xdata = np.linspace(xmin, xmax, 100)
                    ydata = np.zeros(100)
                    plt.plot(xdata, ydata)
                    plt.vlines(tau, 0., P)
                else:
                    print("Multiple pulses list")
                    xmin = self.Pulses[ind].taus[0]   
                    xmax = self.Pulses[ind].taus[-1]    
                    plt.xlim(xmin, xmax)
                    xdata = np.linspace(xmin, xmax, 100)
                    ydata = np.zeros(100)
                    plt.plot(xdata, ydata)
                    for i in range(len(self.Pulses[ind].taus)):
                        tau = self.Pulses[ind].taus[i]
                        P   = self.Pulses[ind].Ps[i] 
                        plt.vlines(tau, 0., P)
            else:
                nPu = len(self.Pulses.taus)
                xmin = self.Pulses.taus[0]   
                xmax = self.Pulses.taus[-1]    
                plt.xlim(xmin, xmax)
                xdata = np.linspace(xmin, xmax, 100)
                ydata = np.zeros(100)
                plt.plot(xdata, ydata, '--')
                taus = self.Pulses.taus
                Ps   = self.Pulses.Ps
                Pmins = np.zeros(nPu)
                plt.vlines(taus, Pmins, Ps, linestyles="solid")
                #for i in range(nPu):
                #    tau = self.Pulses.taus[i]
                #    P   = self.Pulses.Ps[i] 
                #    plt.vlines(tau, 0., P)


        except ind > len(self.Pulses):
            print("Tried to print pulse {}, but there are only {} pulses".format(ind, len(self.Pulses)))






class TwoDTomo(Tomography):
    """
    Class containing methods for 2D tomography common to both impact and Gauss tomography
    """

    def set_istate(self, state):
        """Sets istate, the state to be measured
        Args:
            state : QuTiP Qobj
                New initial state
        """
        if state.dims[0][0] == 2:
            self.istate = state
        else:
            raise Exception("Dimensional mismatch, expected dim = 2, got:", state.dims[0][0])



    @staticmethod
    def get_Stokes_projs(self, rho):
        """Obtain the Stokes paramters, S0, S1, S2, S3 from the density matrix
        Args:
            rho : QuTip Qobj
                The density matrix
        """
        from Utility import Proj
        Stokes = np.zeros(4)
        for i in range(4):
            Op = PauliN(i, self.dim)
            Stokes[i] = 0.5 * Proj(rho, Op)
        return Stokes


    def get_Stokes(self, which='istate') -> None:
        from qutip import ket2dm
        """Obtain the Stokes parameters of the wanted state
        
        Args:
            which : string
                indicates if we are to obtain the Stokes parameters for the 
                initial state, istate, or the recondstructed state, tstate.
        """
        if which == 'istate':
            if self.istate != None:
                if self.istate.type == 'ket':
                    dm = ket2dm(self.istate)
                elif self.istate.type == 'oper':
                    dm = self.istate
                else:
                    raise RuntimeError("Need a ket or density matrix to provide the Stokes' parameters")
                self.iStokes = self.get_2DStokes_projs(self, dm)

        elif which == 'tstate':
            if self.tstate != None:
                if self.tstate.type == 'ket':
                    dm = ket2dm(self.tstate)
                elif self.tstate.type == 'oper':
                    dm = self.tstate
                else:
                    raise RuntimeError("Need a ket or density matrix to provide the Stokes' parameters")
                self.tStokes = self.get_2DStokes_projs(self, dm)

        else:
            raise AttributeError("No valid state was provided")

    def set_Stokes_matrix(self, Projections):
        """Sets the Stokes matrix
        Args: 
            Projections : Projections of $\sigma_3$ on $\sigma_1$, $\sigma_2$ and $\sigma_3$"""
        
        n = 3
        SM = np.zeros((n,n))
        for i in range(n):
            for j in range(n):
                SM[i,j] = Projections[i,j]
        self.StokesMat = Qobj(SM)

    def get_Stokes_matrix(self):
        """Gets the Stokes matrix"""
        from Utility import UBWO, Proj
        if self.StokesMat is not None:
            del self.StokesMat

        Op = 0.5 * (PauliN(0,self.dim) + PauliN(3,self.dim))
        StokesMat = np.zeros((3,3))
        for i in range(3):
            Op_BW = UBWO(self.Us[i].U, Op)
            for j in range(3):
                StokesMat[i,j] = Proj(Op_BW, PauliN(j+1,self.dim))

        self.StokesMat = Qobj(StokesMat)


    def set_S0(self, ProjS0):
        """Sets the projection on the \sigma_0$-matrix
        Args:
            Proj0 : 
                New projection
        """
        self.ProjS0 = ProjS0

    def get_S0(self):
        """Gets the projections of the backwards propagated operator self.Op on the PauliN(0,self.dim) matrix"""
        pass



    def get_tStokes(self, method='straightforward') -> None:
        """Obtaining the tomographically measured Stokes parameters S1, S2, S3
        Args:
            method : string
                Which method to use
        """

        Stokes = np.zeros(4)
        if method == 'straightforward':
            Stokes[0] = 1.
            for i in range(3):
                Stokes[i+1] = 2.*self.pop0s[i] - 1.
            self.tStokes = Stokes
        elif method == 'StokesMat':
            Stokes[0] = 1.
            St = np.zeros(3)
            for i in range(3):
                St[i] = 2.*self.pop0s[i] - 1.
            St = self.StokesMat.inv() * St
            for j in range(3):
                Stokes[j+1] = St[j].real
            self.tStokes = Stokes
        else:
            raise RuntimeError('Method', method, 'not supported yet')



class Tomography2D(ImpactTomography, TwoDTomo):
    """Class handling 2D tomography in the impulse approximation

    Attributes:
    -----------
    Pulses : Pulse parafile 
        Containing the pulse fluencies and time delays
    istate : QuTip Qobj
        Initial state to be measured by tomography
    B : float 
        The rotational constant
    dim : int
        The dimension of the basis
    Projections : 3 by 3 float matrix
        Calculated projections of $\sigma_3$ on \sigma_1, \Sigma_2 and \sigma_3$ after backwards propagation by Pulses
    ProjS0 : float
        Calculated projection on $\sigma_0$ after backwards propagation by Pulses
    Methods : Dictionary
        Information regarding which method was used in order to achieve the Pulses
    Lagrange : Dictionary
        Contains the lagrange multipliers used to optimize the Pulses
    Name : String
        Name of the tomography instance

    Methods:
    --------
    set_pulses(Pulses):
        Sets a pulse with pulse strengths and time delays
    set_EvolutionOperators():
        Sets the three evolution operators used for the tomography scheme

    """

    def __init__(self, Pulses=None, Molpara=None, istate=None, Projections=None, ProjS0=None, Methods=None, Lagrange=None, Name=None):
        """Initialization routine
        Args:
            Pulses (Pulse para): Pulse parafile containing the pulse fluencies and time delays

            B (float): The rotational constant

            dim (int): The dimension of the basis used

            Projections (list of floats): Projections onto the PauliN(3,dim) matrix

            ProjS0 (list of floats): Projections onto the PauliN(0,dim) matrix

            Methods (dict): Which method is used to obtain the tomographic pulses

            Lagrange (dict): The Lagrange multipliers used for optimization

            Name (string): Optional name for the tomography run
        """
        #Add check that we have a 2D initial state
        if istate.dims[0][0] == 2:
            super().__init__(Pulses, Molpara, istate, 2, Name)
            self.iStokes = None
            self.tStokes = None
            self.pop0s   = None

            if Projections is not None:
                self.set_Stokes_matrix(Projections)
            else:
                self.StokesMat = None

            if ProjS0 is not None:
                self.set_S0(ProjS0)
            else:
                self.ProjS0 = None

        else:
            raise Exception("Dimensional mismatch", dim, istate.dims[0][0])


 

    def print_info(self):
        """Print info of the instance of Tomography2D"""
        print(self.Name, "is an instance of", Tomography2D.__name__, ", a", Tomography2D.__doc__)
        print()
        print("Tomography method:", self.Methods)
        print("Basis dimension:", self.dim)
        print("Rotational constant:", self.B, "(au)")
        if self.Pulses is None:
            print("No pulses or delays")
        else:
            for i in range(len(self.Pulses)):
                print("Pulse strengths and time delays:")
                print("P =", self.Pulses[i].P, "(au)", "t = ", self.Pulses[i].t, "(au)")
        if self.StokesMat is None:
            print("No Stokes matrix")
        else:
            print("Stokes matrix:")
            print(self.StokesMat.full())
        if self.Lagrange is None:
            print("No Lagrange multipliers")
        else:
            print("Lagrange multipliers:", self.Lagrange)
        if self.ProjS0 is None:
            print("No S0 projection")
        else:
            print("Projection on S0:", self.ProjS0)



class Tomography2D_BCH(ImpactTomography, TwoDTomo):
    """Class handling 2D tomography in the impulse approximation with corrections from Baker-Campbell-hausdorff formula

    Attributes:
    -----------
    Pulses : Pulse parafile 
        Containing the pulse fluencies and time delays
    istate : QuTip Qobj
        Initial state to be measured by tomography
    B : float 
        The rotational constant
    Projections : 3 by 3 float matrix
        Calculated projections of $\sigma_3$ on \sigma_1, \Sigma_2 and \sigma_3$ after backwards propagation by Pulses
    ProjS0 : float
        Calculated projection on $\sigma_0$ after backwards propagation by Pulses
    Methods : Dictionary
        Information regarding which method was used in order to achieve the Pulses
    Lagrange : Dictionary
        Contains the lagrange multipliers used to optimize the Pulses
    Name : String
        Name of the tomography instance

    Methods:
    --------
    set_pulses(Pulses):
        Sets a pulse with pulse strengths and time delays
    set_EvolutionOperators():
        Sets the three evolution operators used for the tomography scheme

    """
    a_1 = 2. *3.**(-1.) * math.sqrt(5.)**(-1.)
    a_3 = -2. * 21.**(-1.)

    def __init__(self, Pulses=None, Molpara=None, istate=None, Projections=None, ProjS0=None, Methods=None, Lagrange=None, Name=None):
        """Initialization routine
        Args:
            Pulses (Pulse para): Pulse parafile containing the pulse fluencies and time delays

            B (float): The rotational constant

            dim (int): The dimension of the basis used

            Projections (list of floats): Projections onto the PauliN(3,dim) matrix

            ProjS0 (list of floats): Projections onto the PauliN(0,dim) matrix

            Methods (dict): Which method is used to obtain the tomographic pulses

            Lagrange (dict): The Lagrange multipliers used for optimization

            Name (string): Optional name for the tomography run
        """
        super().__init__(Pulses, Molpara, istate, 2, Name)
        try:
            self.Pulses.sigma = Pulses['sigma']
            # Add the BCH-correction to the evolution operators
            self.update_U()
        except KeyError as ke:
            print("KeyError {} in Tomography2D_BCH".format(ke))
        except TypeError as te:
            print("TypeError {} in Tomography2D_BCH".format(te))
        except ValueError as ve:
            print("ValueError {} in Tomography2D_BCH".format(ve))
        except Exception as e:
            print("Exception {} in Tomography2D_BCH".format(e))


    def update_U(self):
        """
        """
        from Utility import U2U1, set_BCH_corr, Uf 
        from qutip import qeye, sigmax, sigmay, sigmaz, Qobj
        from parameters import kappa_from_B_sigma, I0_from_P_sigma
        #kappa = kappa_from_B_sigma(self.Molecule.B, self.Pulses.sigma)
        #kappa = self.Molecule.B * 2. * math.sqrt(2.*math.log(2.)) * self.Pulses.sigma
        #[print("tfact: {}, kappa: {}".format(self.tfact(kappa, P), kappa)) for P in self.Pulses.Ps]
        #thetas = [self.theta_kp(kappa, P) for P in self.Pulses.Ps]
        Id = qeye(2)
        print("Update U")
        for i in range(len(self.Us)):
            P   = self.Pulses.Ps[i]
            tau = self.Pulses.taus[i]
            #U_corr = set_BCH_corr(self.Molecule.B, P, self.Pulses.sigma)
            theta, n0, n1, n2, n3 = set_BCH_corr(self.Molecule.B, P, self.Pulses.sigma)
            pref = cmath.exp(1.j*n0)
            CO   = math.cos(theta) * Id
            SO1  = 1.j * n1 * math.sin(theta) * sigmax()
            SO2  = 1.j * n2 * math.sin(theta) * sigmay()
            SO3  = 1.j * n3 * math.sin(theta) * sigmaz()
            Upulse = pref * (CO + SO1 + SO2 + SO3)
            fwhm = 2.*math.sqrt(2.*math.log(2.)) * self.Pulses.sigma
            Ufree  = Uf(self.Molecule.B, tau-fwhm, 2)
            self.Us[i].U = 0. * U2U1(Upulse, Ufree)





    def tfact(self, kappa, P):
        """
        """
        f1 = 3.*kappa * (3.*kappa/self.a_1 - 2.*P*self.a_3)
        f2 = P**2.*(self.a_1**2. + self.a_3**2.) +9.
        return math.sqrt(f1 + f2)


    def theta_kp(self, kappa, P):
        fact = self.tfact(kappa, P)
        return self.a_1 * kappa * P * fact

    def theta_kp_deg(self, kappa, P):
        return self.theta_kp(kappa, P) * 180. / math.pi

    def n_1_kp(self, kappa, P):
        f1 = math.sqrt(3.*kappa*(3.*kappa/self.a_1 - 2.*P*self.a_3) + (P*self.a_3)**2.)
        f2 = self.tfact(kappa ,P)
        return f1 * f2**(-1.)

    def n_2_kp(self, kappa, P):
        f1 = 3.
        f2 = self.tfact(kappa, P)
        return f1 * f2**(-1.)

    def n_3_kp(self, kappa, P):
        f1 = P * self.a_1
        f2 = self.tfact(kappa, P)
        return f1 * f2**(-1.)



class Tomography4D_BCH(ImpactTomography):
    """Class handling 2D tomography in the impulse approximation with corrections from Baker-Campbell-hausdorff formula

    Attributes:
    -----------
    Pulses : Pulse parafile 
        Containing the pulse fluencies and time delays
    istate : QuTip Qobj
        Initial state to be measured by tomography
    B : float 
        The rotational constant
    Projections : 3 by 3 float matrix
        Calculated projections of $\sigma_3$ on \sigma_1, \Sigma_2 and \sigma_3$ after backwards propagation by Pulses
    ProjS0 : float
        Calculated projection on $\sigma_0$ after backwards propagation by Pulses
    Methods : Dictionary
        Information regarding which method was used in order to achieve the Pulses
    Lagrange : Dictionary
        Contains the lagrange multipliers used to optimize the Pulses
    Name : String
        Name of the tomography instance

    Methods:
    --------
    set_pulses(Pulses):
        Sets a pulse with pulse strengths and time delays
    set_EvolutionOperators():
        Sets the three evolution operators used for the tomography scheme

    """
    #a_1 = 2. *3.**(-1.) * math.sqrt(5.)**(-1.)
    #a_3 = -2. * 21.**(-1.)

    def __init__(self, Pulses=None, Molpara=None, istate=None, Projections=None, ProjS0=None, Methods=None, Lagrange=None, Name=None):
        """Initialization routine
        Args:
            Pulses (Pulse para): Pulse parafile containing the pulse fluencies and time delays

            B (float): The rotational constant

            dim (int): The dimension of the basis used

            Projections (list of floats): Projections onto the PauliN(3,dim) matrix

            ProjS0 (list of floats): Projections onto the PauliN(0,dim) matrix

            Methods (dict): Which method is used to obtain the tomographic pulses

            Lagrange (dict): The Lagrange multipliers used for optimization

            Name (string): Optional name for the tomography run
        """
        super().__init__(Pulses, Molpara, istate, 4, Name)
        try:
            self.Pulses.sigma = Pulses['sigma']
            # Add the BCH-correction to the evolution operators
            self.update_U()
        except KeyError as ke:
            print("KeyError {} in Tomography2D_BCH".format(ke))
        except TypeError as te:
            print("TypeError {} in Tomography2D_BCH".format(te))
        except ValueError as ve:
            print("ValueError {} in Tomography2D_BCH".format(ve))
        except Exception as e:
            print("Exception {} in Tomography2D_BCH".format(e))



    def update_U(self):
        """
        """
        from Utility import U2U1, set_BCH_corr4D 
        from qutip import qeye, sigmax, sigmay, sigmaz, Qobj
        from parameters import kappa_from_B_sigma, I0_from_P_sigma
        
        for i in range(len(self.Us)):
            P = self.Pulses.Ps[i]
            U_corr = set_BCH_corr4D(self.Molecule.B, P, self.Pulses.sigma)
            # Need to multiply the data togenter because the tensor product causes dimensional problems   
            U_data = self.Us[i].U * U_corr.data# * self.Us[i].U.data
            self.Us[i].U = Qobj(U_data) #U2U1(U_corr, self.Us[i].U)





class Tomography3D(ImpactTomography):
    """Class handling 3D tomography in the impulse approximation

    Attributes:
    -----------
    Pulses : Pulses parafile
        Containing the pulse fluencies and time delays
    istate : QuTip Qobj
        Initial state to be measured by tomography
    B : float 
        The rotational constant
    dim : int
        The dimension of the basis
    Projections : 3 by 3 float matrix
        Calculated projections of $\sigma_3$ on \sigma_1, \Sigma_2 and \sigma_3$ after backwards propagation by Pulses
    ProjS0 : float
        Calculated projection on $\sigma_0$ after backwards propagation by Pulses
    Methods : Dictionary
        Information regarding which method was used in order to achieve the Pulses
    Lagrange : Dictionary
        Contains the lagrange multipliers used to optimize the Pulses
    Name : String
        Name of the tomography instance

    Methods:
    --------
    set_pulses(Pulses):
        Sets a pulse with pulse strengths and time delays
    set_EvolutionOperators():
        Sets the three evolution operators used for the tomography scheme

    """

    def __init__(self, Pulses=None, istate=None, B=1., dim=4, Projections=None, ProjS0=None, Methods=None, Lagrange=None, Name=None):
        """Initialization routine
        Args:
            Pulses (Pulsepara): Pulses parafile containing the pulse fluencies and time delays

            B (float): The rotational constant

            dim (int): The dimension of the basis used

            Projections (list of floats): Projections onto the PauliN(3,dim) matrix

            ProjS0 (list of floats): Projections onto the PauliN(0,dim) matrix

            Methods (dict): Which method is used to obtain the tomographic pulses

            Lagrange (dict): The Lagrange multipliers used for optimization

            Name (string): Optional name for the tomography run
        """
        #Add check that we have a 3D initial state
        if (dim >= 3) and (istate.dims[0][0] == 3):
            super().__init__(Pulses, istate, B, dim, Name)
            self.iStokes = None
            self.tStokes = None
            self.pop0s   = None

            if Projections is not None:
                self.set_Stokes_matrix(Projections)
            else:
                self.StokesMat = None

            if ProjS0 is not None:
                self.set_S0(ProjS0)
            else:
                self.ProjS0 = None

        else:
            raise Exception("Dimensional mismatch", dim, istate.dims[0][0])


    def set_istate(self, state):
        """Sets istate, the state to be measured
        Args:
            state : QuTiP Qobj
                New initial state
        """
        if state.dims[0] == 3:
            self.istate = state
        else:
            raise ValueError("Dimensional mismatch, expected dim = 3, got:", state.dims[0])


    @staticmethod
    def get_3DStokes_projs(self, rho):
        """Obtain the Stokes paramters, S0, S1, S2, S3 from the density matrix
        Args:
            rho : QuTip Qobj
                The density matrix
        """
        Stokes = np.zeros(9)
        for i in range(9):
            Op = GMN(i, self.dim)
            Stokes[i] = Proj(rho, Op)
        return Stokes


class GaussTomography(Tomography):
    """
    Tomography using real Gaussian pulses
    """

    def __init__(self, Pulsepara=None, Molpara=None, istate=None, dim=3, Methods=None, Lagrange=None, Name=None):
        """
        Constructor
        """
        super().__init__(Pulsepara, Molpara, istate, dim, Methods, Lagrange, Name, Ptype="Gauss")
        self.Da   = self.Molecule.Da #Molpara['Da'] Remove this later
        #self.set_pulses(Pulsepara)
        self.set_Hams()
        #self.H0   = None
        #self.Hint = None


    def set_pulses(self, Pulses, Name=None) -> None:
        """Sets the pulses for the tomography instance
        
        Args: 
            Pulses : Pulsepara containing the pulses strengths and time delays

        """
        import Utility as Ut
        nP  = len(Pulses['I0s'])
        #self.Pulses  = []
        #for i in range(nP):
        #    Para = {'I0s': [Pulses['I0s'][i]], 'taus': [Pulses['taus'][i]], 'sigma': Pulses['sigma'], 'tmax': Pulses['tmax'], 'tmin': Pulses['tmin'], 'dt': Pulses['dt']}
            #Pulse = Ut.Pulses(Pulses[i].P, Pulses[i].t)
            #Pulse = Ut.GaussPulses(Para,  Name)
            #self.Pulses.append(Pulse)
        self.Pulses = Ut.GaussPulses(Pulses)


    def set_H0(self):
        """
        """
        from Utility import H0
        self.H0   = H0(self.B, self.dim)

    def set_Hint(self):
        """
        """
        from Utility import H1_I0
        self.Hint = H1_I0(self.Da, self.dim)


    def set_Hams(self):
        """
        """
        self.set_H0()
        self.set_Hint()


    def prop_FW(self, args, funct, options=None):
        """
        Propagate forward under the influence of the Gaussian pulse
        """
        from qutip import sesolve
        from Utility import Gauss_me, double_Gauss_me, N_Gauss_me

        res = sesolve([self.H0, [self.Hint, Gauss_me]], self.istate, np.arange(self.Pulses.tmin, self.Pulses.tmax, self.Pulses.dt), args=args, options=options)


        return res.states[-1]

    def prop_BW_state(self, args, istate, funct, options=None):
        """
        Propagate backwards under the influence of the Gaussian pulse
        """
        from qutip import sesolve
        from Utility import Gauss_me, double_Gauss_me, N_Gauss_me
        from qutip import sigmaz
        

        res = sesolve([self.H0, [self.Hint, funct]], istate, np.arange(self.Pulses.tmax, self.Pulses.tmin, -self.Pulses.dt), args=args, options=options)

        return res.states[-1]



    def prop_BW(self, args, funct, options=None):
        """
        Propagate the extended \sigma3 operator backwards under the influence of the Gaussian pulse
        """
        from qutip import sesolve, basis
        from Utility import Gauss_me, double_Gauss_me, N_Gauss_me
        from qutip import ket2dm
        

        res1 = sesolve([self.H0, [self.Hint, funct]], basis(self.dim, 0), np.arange(self.Pulses.tmax, self.Pulses.tmin, -self.Pulses.dt), args=args, options=options)

        res2 = sesolve([self.H0, [self.Hint, funct]], basis(self.dim, 1), np.arange(self.Pulses.tmax, self.Pulses.tmin, -self.Pulses.dt), args=args, options=options)

        #sig3_BW = res1.states[-1] * res1.states[-1].dag() - res2.states[-1] * res2.states[-1].dag()
        sig3_BW = ket2dm(res1.states[-1]) - ket2dm(res2.states[-1])

        return sig3_BW

    def pop_runs(self, options=None):
        """
        """
        from qutip import Options
        if options is None:
            options = Options()
        pops = np.zeros(self.Pulses.nP)
        for i in self.Pulses.nP:
            args = {'I0s': self.Pulses.I0s[i], 'taus': self.Pulses.taus[i], 'sigma': self.Pulses.sigma}
            pop = self.prop_FW(args, options)
            pops[i] = pop

        self.pops = pops


    def tomo_runs(self, options):
        """
        """
        if options is None:
            options = Options()
        Ops = []
        for i in range(self.Pulses.nP):
            args = {'I0s': self.Pulses.I0s[i], 'taus': self.Pulses.taus[i], 'sigma': self.Pulses.sigma}
            Op = self.prop_BW(args, options)
            Ops.append(Op)

        self.Ops = Ops


    def optimize(self, Op, variables:dict, params:dict, method:dict, target:dict, lagrange:float):
        """
        We want to optimize pulse parameters to implement dynamics of a given measurment operator.

        Args:
        -----
        
        """
        from Opt import GaussOptimizer
        Opt = GaussOptimizer(Op, self.dim, variables, params, method, target, lagrange)
        Opt.optimization()

#    def print_info(self, full_mol_info=False):
#        """
#        """
#        if self.Name is None:
#            print("Printing info on nameless {}Tomography instance".format(self.Ptype))
#        else:
#            print("Printing info for {}Tomography instance called".format(self.Ptype, self.Name))
#        print()
#        print("Molecule info:")
#        print("--------------")
#        self.Molecule.print_info(full_mol_info)
#        print()
#        print("Pulse info:")
#        print("-----------")
#        #for i in range(len(self.Pulses)):
#        #    self.Pulses[i].print_pulse_info()
#        self.Pulses.print_info()



class GaussTomography2D(GaussTomography, TwoDTomo):
    """
    A class for performing Tomography with a Gaussian pulse in a 2D basis.
    """
    def __init__(self, Pulsepara=None, Molpara=None, istate=None, Methods=None, Lagrange=None, Name=None):
        """
        """
        try:
            super().__init__(Pulsepara=Pulsepara, Molpara=Molpara, istate=istate, dim=2, Methods=Methods, Lagrange=Lagrange, Name=Name)
            self.set_istate(istate)
        except KeyError as ke:
            print("KeyError: {} occurred in GaussTomography2D".format(ke))
        except TypeError as te:
            print("TypeError: {} occurred in GaussTomography2D".format(te))
        except ValueError as ve:
            print("ValueError: {} occurred in GaussTomography2D".format(ve))
        except Exception as e:
            print("Exception: {} occurred in GaussTomography2D".format(e))
        


    def run_tomography():
        """
        Runs a tomography scheme to retrieve the initial state istate
        """
    
        # We begin by setting the hamiltonians (if they don't already exist)
        if self.H0 is None:
            self.set_H0()

        # Noe we set the pulses
        self.set_pulses()

        # Then we run the populaton runs 
        self.pop_runs()

        # Finally we do the backwards tomography runs
        self.tomo_runs()



def Domega(P):
    Do = -P * sqdisc
    return Do

def Stokesi(i, P):

    if i == 0:
        return 1.
    elif i == 1:
        DO = Domega(P)
        return 2. * aplus * c2b * (aplus**2. - c2b**2.) * (1. - math.cos(DO)) / nplus**4. 
    elif i == 2:
        DO = Domega(P)
        return 2. * aplus * c2b * (aplus**2. + c2b**2.) * math.sin(DO) / nplus**4.
    elif i == 3:
        DO = Domega(P)
        return ( (aplus**2. - c2b**2.)**2. + 4. * (aplus * c2b)**2. * math.cos(DO)) / nplus**4.
    else:
        raise ValueError("Index", i, "out of range for Stokes: 0,1,2,3")

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
    T2 = math.pi/(12.*B) + TS1(B, n)
    
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
            Po = sigmaz()
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


