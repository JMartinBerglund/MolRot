"""
This is the Tomography module
"""

#!/usr/bin/env python3
from textwrap import dedent
import os
import sys
sys.path.append("/home/martin/Work/GitHub/MolRot/Utility/")
sys.path.append("/home/martin/Work/GitHub/MolRot/Tomography/")
import numpy as np
import math
import cmath
from qutip import *
from parameters import eps0, c, fs2au
from Tomomod import *

class Tomography():
    """Class handling tomography in the impulse approximation

    Attributes:
    -----------
    Pulses : Pulses object 
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


    def __init__(self, Pulses=None, istate=None, B=1., dim=3, Methods=None, Lagrange=None, Name=None) -> None:
        """Initialization routine
        
        Args:
            Pulses: Utility Pulses object
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
        self.dim      = dim
        self.istate   = istate
        self.tstate   = None
        self.B        = B
        self.Methods  = Methods
        self.Lagrange = Lagrange
        self.Name     = Name
        
        if Pulses is not None:
            self.set_pulses(Pulses)
            self.set_EvolutionOperators()
        else:
            self.Pulses = None
            self.Us = None

 
    def set_pulses(self, Pulses) -> None:
        """Sets the pulses for the tomography instance
        
        Args: 
            Pulses : Array of pulses instance containing the pulses strengths and time delays
        """
        import Utility as Ut
        nP  = len(Pulses)
        self.Pulses  = []
        for i in range(nP):
            Pulse = Ut.Pulses(Pulses[i].P, Pulses[i].t)
            self.Pulses.append(Pulse)

    def set_EvolutionOperators(self) -> None:
        """Sets the three evolution operators for the tomography scheme"""
        from Utility import EvolutionOperators
        nP  = len(self.Pulses)
        print("Number of pulses:", nP)
        self.Us  = []
        for i in range(nP):
            self.Us.append(EvolutionOperators(self.Pulses[i], self.B, dim=self.dim))


class Tomography2D(Tomography):
    """Class handling 2D tomography in the impulse approximation

    Attributes:
    -----------
    Pulses : Pulses object 
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

    def __init__(self, Pulses=None, istate=None, B=1., dim=3, Projections=None, ProjS0=None, Methods=None, Lagrange=None, Name=None) -> None:
        """Initialization routine
        Args:
            Pulses (Pulses object): Pulses object containing the pulse fluencies and time delays

            B (float): The rotational constant

            dim (int): The dimension of the basis used

            Projections (list of floats): Projections onto the PauliN(3,dim) matrix

            ProjS0 (list of floats): Projections onto the PauliN(0,dim) matrix

            Methods (dict): Which method is used to obtain the tomographic pulses

            Lagrange (dict): The Lagrange multipliers used for optimization

            Name (string): Optional name for the tomography run
        """
        #Add check that we have a 2D initial state
        if (dim >= 2) and (istate.dims[0][0] == 2):
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
        if state.dims[0] == 2:
            self.istate = state
        else:
            raise Exception("Dimensional mismatch, expected dim = 2, got:", state.dims[0])

    @staticmethod
    def get_Stokes_projs(self, rho):
        """Obtain the Stokes paramters, S0, S1, S2, S3 from the density matrix
        Args:
            rho : QuTip Qobj
                The density matrix
        """

        Stokes = np.zeros(4)
        for i in range(4):
            Op = PauliN(i, self.dim)
            Stokes[i] = Proj(rho, Op)

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
                    print("Need a ket or density matrix to provide the Stokes' parameters")
                    raise Exception
                self.iStokes = self.get_Stokes_projs(self, dm)

        elif which == 'tstate':
            if self.tstate != None:
                if self.tstate.type == 'ket':
                    dm = ket2dm(self.tstate)
                elif self.tstate.type == 'oper':
                    dm = self.tstate
                else:
                    print("Need a ket or density matrix to provide the Stokes' parameters")
                    raise Exception
                self.tStokes = self.get_Stokes_projs(self, dm)

        else:
            raise Exception

    def set_Stokes_matrix(self, Projections):
        """Sets the Stokes matrix
        Args: 
            Projections : Projections of $\sigma_3$ on $\sigma_1$, $\sigma_2$ and $\sigma_3$
        """
        
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
            print('Method', method, 'not supported yet')
            raise Exception

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


