"""Module containing interferometry classes and methods in the impulse approximation and without the approximation"""

#!/usr/bin/env python3
from textwrap import dedent
import os
import sys
sys.path.append("/home/martin/Work/quantumLibWorkP3")
import numpy as np
import math
import cmath
from qutip import *
import parameters as pm
from Tensor import RotationTensor


"""
Some useful interferometry methods
"""

def Efield_interferometry_impact_2level(Brot, eps0, D, P, tau):
    """
    Runs an analytical interferometry aimed at measuring the local E-field based in the impact approximation in 2-level (2-even and 2 odd levels)
    """
    from Utility import MatrixElementsPulseImpact2D as ME
    from Utility import EigenField2D as EF2D
    #print("WARNING NOT PROPERLY TESTED YET!!!")

    # Obtaining the nececcary matrix elements
    A, B,_ = ME(P)
    E,_,_  = ME(P, odd=True)

    # Interaction strength and eigenvalues of the dipole interaction
    l_m, l_p, V = EF2D(Brot, D, eps0)
    Dl  = 2. * math.sqrt(Brot**2. + V**2.)

    # Free egienvalues of the states j=1,2
    omega1 = 2.* Brot
    omega2 = 6.* Brot

    # Prefactor
    fact = math.sqrt(l_p**2. + V**2.)**(-1.)

    # The Non-oscillating term X0
    f00 = A * (A*l_p**2. + E*V**2.)
    f01 = A * V**2. * (A - E)
    X0  = l_p**2. * (fact**6. * (f00 * np.conj(f00) + f01 * np.conj(f01))  + fact**2. *(B*np.conj(B))**2. ) 

    # Now for the oscilating terms
    lt = len(tau)
    c_sq = np.zeros(lt)
    for i in range(lt):
        # The term oscillating at frequency Dl
        e1  = cmath.exp(1.j*Dl*tau[i])
        f1  = f00 * np.conj(f01)
        f11 = f1 * e1
        X1  = 2. * l_p**2. * fact**6. * f11.real

        # The term oscillating at frequency omega2 - l_m
        e2  = cmath.exp(1.j * (omega2 - l_m) * tau[i])
        f2  = f00 * np.conj(B**2.)
        f22 = f2 * e2
        X2  = 2. * l_p**2. * fact**4. * f22.real

        # The term oscillating at frequency omega2 - l_p
        e3  = cmath.exp(1.j * (omega2 - l_p) * tau[i])
        f3  = f01 * np.conj(B**2.)
        f33 = f3 * e3
        X3  = 2. * l_p**2. * fact**4. * f33.real

        c_sq[i] = X0 + X1 + X2 + X3

    return c_sq

def get_dtau(omega_max: float) -> float:
    """
    Obtain the minimum time spacing (Delta tau) needed to resolve the maximum frequency given. (From the uncertainty principle)

    Args:
    -----
        omega_max: float
            The maximum frequency we want to resolve

    Returns:
    --------
        dtau_min: float
            The minimum time spacing needed
    """

    return 0.5 * omega_max**(-1.)

def get_taumax(domega:float) -> float:
    """
    """
    return 0.5 * domega**(-1.)

def get_ntau(taumax:float, dtau:float) -> float:
    """
    """
    return int((taumax + dtau) / dtau)

def taus_from_mol(Molecule, omax, do):
    Mol       = pm.set_MolParams(Molecule)
    omega_max = omax * Mol['B']
    domega    = do * Mol['B']
    dtau      = get_dtau(omega_max)
    taumax    = get_taumax(domega)
    ntau      = get_ntau(taumax, dtau)
    
    return taumax, dtau, ntau





class Interferometry():
    """
    The main interferometry class
    """

    #def __init__(self, istate, mstate=0, B=1., dim=3, a=0., Name=None, odd=False, mrep=False) -> None:
    def __init__(self, Molecule, istate, mstate=0, dim=3, Name=None, odd=False, mrep=False, custom_para=None) -> None:
        """
        Initializes an instance of the interferometer
            
            Args:
            -----
                Molecule: dict|string
                    Contains the molecular paramters B and optionally a
                    B: float
                        The rotational constant
                    a: float
                        Centrifugal distortional constant
                    Mol: string
                        Which molecule to use. If only string, can be MgH+, CaH+, CaD+, Custom

                istate: Qobj (ket or density matrix)
                    Initial state for the interferometer. 
                
                mstate : int
                    Which final state to measure the interferogram on

                dim: int
                    The dimansion of the basis used

                Name: str
                    Name of the interferomoeter intance

                mrep: Boolean
                    If true, include non-zero m-states in the representation

                kwargs: dictionary
                    If Molecule is Custom, dictioinary for initialization

        """
        from parameters import Molecule as M
        try:
            if not isinstance(Molecule, str):
                self.check_para(Molecule)
            self.molecule = M(Molecule, custom_para)
            self.mstate   = mstate
            self.dim      = dim
            if Name is None:
                self.Name = None
            else:
                self.Name = Name
            self.U        = None
            self.inter    = None
            self.tau      = None
            self.spectrum = None
            self.omega    = None
            self.odd      = odd
            self.mrep     = mrep
            self.initialize_istate(istate)

        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except ValueError as ve:
            print("ValueError: {}".format(ve))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except Exception as e:
            print("Exception: {}".format(e))

    def check_para(self, para:dict|str):
        """
        Checking the paramter dictionary for consistency.
        """
        if isinstance(para, dict):

            if 'B' in para:
                if isinstance('B', float):
                    if 'B' < 0.:
                        raise ValueError("The provided rotational constant is negative: {}".format(para['B']))
                else:
                    raise TypeError("Expected a real rotational constant, the provided one is: {}".format(type(para['B'])))
            else:
                raise KeyError("No rotational constant provided in the dictionary")

            if 'a' in para:
                if isinstance('a', float):
                    if 'a' > 0.:
                        raise ValueError("The provided rotational distorsional constant is positive: {}".format(para['a']))
                else:
                    raise TypeError("Expected a real rotational distorsional constant, the provided one is: {}".format(type(para['a'])))
            else:
                print("No distorsional constant provided, setting to 0")
                para.update({'a': 0.})
        

    def initialize_istate(self, istate):
        """
        Initializes the given state
        """
        #from parameters import set_state
        #self.istate = set_state(self.istate, self.dim, self.mrep)
        if isinstance(istate, int):
            if self.mrep:
                if istate < (self.dim + 1)**2:
                    self.istate = qutip.basis((self.dim + 1)**2,istate)
                else:
                    raise ValueError("Can't set to a state larger than the basis size")
            else:
                if istate < self.dim:
                    self.istate = qutip.basis(self.dim, istate)
                else:
                    raise Exception("Can't set to a state larger than the basis size")

        elif hasattr(istate, 'type'):
            if istate.type == 'ket':
                self.istate = istate
            else:
                raise TypeError("A Qobject must be of ket type, recieved: {}".format(istate.type))
        else:
            print("Warning, no initial state was provided, setting to ground state!")
            if self.mrep:
                self.istate = qutip.basis((self.dim + 1)**2,0)
            else:
                self.istate = qutip.basis(self.dim, 0)


    def print_info(self) -> None:
        """
        Printing info about the interferometry object.
        """
        print("Info about interferometry:")
        print("--------------------------")
        print("")
        print("Molecular system:")
        print("--------------------------")
        self.molecule.print_info()

    def run_interferometry(self) -> None:
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        pass

    def get_visibility(self):
        """
        """
        if self.inter is None:
            raise TypeError("The interferogram is none.")
        else:
            Vplus  = max(self.inter)
            Vminus = min(self.inter)
            return (Vplus - Vminus) / (Vplus + Vminus)
        
    def run_2levelImpactEfield_interferometry(self, eps0, D, P, tau) -> None:
        """
        Runs an analytical interferometry aimed at measuring the local E-field based in the impact approximation in 2-level (2-even and 2 odd levels)
        """
        from Utility import MatrixElementsPulseImpact2D as ME
        from Utility import EigenField2D as EF2D
        #print("WARNING NOT PROPERLY TESTED YET!!!")

        # Obtaining the nececcary matrix elements
        A, B,_ = ME(P)
        E,_,_  = ME(P, odd=True)

        # Interaction strength and eigenvalues of the dipole interaction
        l_m, l_p, V = EF2D(self.B, D, eps0)
        Dl  = 2. * math.sqrt(self.B**2. + V**2.)

        # Free egienvalues of the states j=1,2
        omega1 = 2.* self.B
        omega2 = 6.* self.B

        # Prefactor
        fact = math.sqrt(l_p**2. + V**2.)**(-1.)

        # The Non-oscillating term X0
        f00 = A * (A*l_p**2. + E*V**2.)
        f01 = A * V**2. * (A - E)
        X0  = l_p**2. * (fact**6. * (f00 * np.conj(f00) + f01 * np.conj(f01))  + fact**2. *(B*np.conj(B))**2. ) 

        # Now for the oscilating terms
        lt = len(tau)
        c_sq = np.zeros(lt)
        for i in range(lt):
            # The term oscillating at frequency Dl
            e1  = cmath.exp(1.j*Dl*tau[i])
            f1  = f00 * np.conj(f01)
            f11 = f1 * e1
            X1  = 2. * l_p**2. * fact**6. * f11.real

            # The term oscillating at frequency omega2 - l_m
            e2  = cmath.exp(1.j * (omega2 - l_m) * tau[i])
            f2  = f00 * np.conj(B**2.)
            f22 = f2 * e2
            X2  = 2. * l_p**2. * fact**4. * f22.real

            # The term oscillating at frequency omega2 - l_p
            e3  = cmath.exp(1.j * (omega2 - l_p) * tau[i])
            f3  = f01 * np.conj(B**2.)
            f33 = f3 * e3
            X3  = 2. * l_p**2. * fact**4. * f33.real

            c_sq[i] = X0 + X1 + X2 + X3

        return c_sq

    @staticmethod
    def opt_2levelImpactEfield_interferometry(Brot, eps0, D, P, tau) -> None:
        """
        Runs an analytical interferometry aimed at measuring the local E-field based in the impact approximation in 2-level (2-even and 2 odd levels)
        """
        from Utility import MatrixElementsPulseImpact2D as ME
        from Utility import EigenField2D as EF2D
        print("WARNING NOT PROPERLY TESTED YET!!!")

        # Obtaining the nececcary matrix elements
        A, B,_ = ME(P)
        E,_,_  = ME(P, odd=True)

        # Interaction strength and eigenvalues of the dipole interaction
        l_m, l_p, V = EF2D(Brot, D, eps0)
        Dl  = 2. * math.sqrt(Brot**2. + V**2.)

        # Free egienvalues of the states j=1,2
        omega1 = 2.* Brot
        omega2 = 6.* Brot

        # Prefactor
        fact = math.sqrt(l_p**2. + V**2.)**(-1.)

        # The Non-oscillating term X0
        f00 = A * (A*l_p**2. + E*V**2.)
        f01 = A * V**2. * (A - E)
        X0  = l_p**2. * (fact**6. * (f00 * np.conj(f00) + f01 * np.conj(f01))  + fact**2. *(B*np.conj(B))**2. ) 

        # Now for the oscilating terms
        lt = len(tau)
        c_sq = np.zeros(lt)
        for i in range(lt):
            # The term oscillating at frequency Dl
            e1  = cmath.exp(1.j*Dl*tau[i])
            f1  = f00 * np.conj(f01)
            f11 = f1 * e1
            X1  = 2. * l_p**2. * fact**6. * f11.real

            # The term oscillating at frequency omega2 - l_m
            e2  = cmath.exp(1.j * (omega2 - l_m) * tau[i])
            f2  = f00 * np.conj(B**2.)
            f22 = f2 * e2
            X2  = 2. * l_p**2. * fact**4. * f22.real

            # The term oscillating at frequency omega2 - l_p
            e3  = cmath.exp(1.j * (omega2 - l_p) * tau[i])
            f3  = f01 * np.conj(B**2.)
            f33 = f3 * e3
            X3  = 2. * l_p**2. * fact**4. * f33.real

            c_sq[i] = X0 + X1 + X2 + X3

        return c_sq


    def get_dtau(omega_max: float) -> float:
        """
        Obtain the minimum time spacing (Delta tau) needed to resolve the maximum frequency given. (From the uncertainty principle)

        Args:
        -----
            omega_max: float
                The maximum frequency we want to resolve

        Returns:
        --------
            dtau_min: float
                The minimum time spacing needed
        """

        return 0.5 * omega_max**(-1.)

    def get_taumax(domega:float) -> float:
        """
        """
        return 0.5 * domega**(-1.)

    def get_ntau(taumax:float, dtau:float) -> float:
        """
        """
        return int((taumax + dtau) / dtau)

    def taus_from_mol(Molecule, omax, do):
        Mol       = pm.set_MolParams(Molecule)
        omega_max = omax * Mol['B']
        domega    = do * Mol['B']
        dtau      = get_dtau(omega_max)
        taumax    = get_taumax(domega)
        ntau      = get_ntau(taumax, dtau)
    
        return taumax, dtau, ntau




    def get_spectrum(self, shift=False) -> None:
        """
        Obtains the Fourier spectrum of the interferogram
        """
        from scipy.fft import fft, fftfreq, fftshift
        # Check that the interferogram and delay arrays exist
        # and take the Fourier transform if they do
        if (self.Pulses.taus is not None) and (self.inter is not None):
            ntau = len(self.Pulses.taus)
            dtau = self.Pulses.taus[1] - self.Pulses.taus[0]
            Fourier = fft(self.inter)
            if shift:
                self.spectrum = fftshift(np.abs(Fourier)) / ntau
                self.omega = fftshift(fftfreq(ntau, dtau))
            else:
                self.spectrum = np.abs(Fourier[0:ntau//2]) / ntau
                self.omega = fftfreq(ntau, dtau)[0:ntau//2]

        else:
            print("Tried to take the Fourier transform of interferogram. Can't do that if either the interferogram or time grid don' t exist.")
            raise Exception()
        
    def plot_interferogram(self, title=None, scaleTrot=False, xmin=None, xmax=None, ymin=None, ymax=None) -> None:
        import matplotlib.pyplot as plt
        if (self.Pulses.taus is not None) and (self.inter is not None):
            # Prepare plot
            if scaleTrot:
                Trot = 2.*math.pi * self.molecule.B**(-1.)
                plt.plot(self.Pulses.taus / Trot, self.inter)
            else:
                plt.plot(self.Pulses.taus, self.inter)

            # Set axis labels
            if not scaleTrot:
                plt.xlabel('$\\tau$ (atomic units)')
            else:
                plt.xlabel('$\\tau$ ($T_{rot}^{-1}$)')
            plt.ylabel('Interferogram')

            # Set the title
            if title is not None:
                plt.title(title)
            
            # Scaling
            if xmin is not None:
                plt.xlim(left=xmin)
            if xmax is not None:
                plt.xlim(right=xmax)
            if ymin is not None:
                plt.ylim(bottom=ymin)
            if ymax is not None:
                plt.ylim(top=ymax)

            # Plot the graph
            plt.plot()
        else:
            print("Tried to plot the Fourier spectrum. Can't do that if either the spectrum or omega grid don' t exist.")
            raise Exception()

    def plot_spectrum(self, title=None, scaleTrot=False, setTicks=False, xmin=None, xmax=None, ymin=None, ymax=None) -> None:
        import matplotlib.pyplot as plt
        if (self.omega is not None) and (self.spectrum is not None):
            # Prepare plot
            if scaleTrot is True:
                # Change to get_Trot()
                Trot = 2.*math.pi*(self.molecule.B)**(-1.)
                plt.plot(self.omega * Trot, self.spectrum)
            else:
                plt.plot(self.omega, self.spectrum)


            # Set axis labels 
            if scaleTrot:
                plt.xlabel('$\omega$ ($T_{rot}$)' )
                if setTicks:
                    plt.xticks(np.arange(self.omega[0]*Trot, self.omega[-1]*Trot, step=2))
            else:
                plt.xlabel('$\omega$ (au)')

                    
            plt.ylabel('Spectrum')

            # Set the title
            if title is not None:
                plt.title(title)

            # Scaling
            if xmin is not None:
                plt.xlim(left=xmin)
            if xmax is not None:
                plt.xlim(right=xmax)
            if ymin is not None:
                plt.ylim(bottom=ymin)
            if ymax is not None:
                plt.ylim(top=ymax)

            # Plot the graph
            plt.plot()
            plt.show()
        else:
            print("Tried to plot the Fourier spectrum. Can't do that if either the spectrum or omega grid don' t exist.")
            raise Exception()

    @classmethod
    def get_omega0(self, omega, spec, p0=None):
        """
        Perfroming the curve fit to get the energy difference between the (0,0) and (1,0) states
        """
        from scipy.optimize import curve_fit as cf
        import Utility as Ut

        if p0 is None:
            [param, covar] = cf(Ut.Lorentz, omega, spec)
        else:
            [param, covar] = cf(Ut.Lorentz, omega, spec, p0)

        return param

    
    def find_omega0(self, omin=None, omax=None, p0=None, om1 = 200.):
        """
        Obtain the measured energy difference between the (0,0) and (1,0) states from the spectrum by curve fit
        """
        if omin is None:
            omin = 1
        if omax is None:
            omax = len(self.omega) - 1
      
        #print(self.get_Trot())
        Trot = 2.*math.pi*self.molecule['B']**(-1.)
        # Select the interval over wich to measure
        om = self.omega[omin:omax] * Trot
        spec = self.spectrum[omin:omax]
        oml = len(spec)
        omega = np.arange(om[0], om[-1], (om[-1] - om[0]) / oml)
        omega_inter = np.linspace(om[0], om[-1], int(om1))
        print(len(spec), len(omega), omega_inter[0], omega_inter[-1])
        # Interpolate the spectrum
        spectrum = np.interp(omega_inter, omega, spec)
        #Finally get the energy difference from curve fit
        omega0 = self.get_omega0(omega_inter, spectrum, p0)
        
        return omega0

    @classmethod
    def get_Trot(self):
        """
        Obtain the rotational period from the rotational constant
        """

        return 2.*math.pi * self.molecule['B']**(-1.)

    def eps2D(self, Deltalambda):
        """
        Obtain the 2D approximation of the electric field 
        """
        #Dl = Deltalambda / self.get_Trot()
        try:
            return self.molecule['D']**(-1) * math.sqrt(3. * (0.25*Deltalambda**2. - self.molecule['B']**2.))
        except KeyError as ke:
            print("KeyError: {} not in molecule".format(ke))
        except ValueError as ve:
            print("ValueError: The square root is not defined for the imput values, resulted in {} under the square root".format(ve))

    @staticmethod
    def find_Efield(x, B, D, dim, deltalambda):
        import Utility as Ut
        eps = x[0]
        print(eps)

        H0 = Ut.H0(B, dim, full=True) + Ut.H_dip(eps, D, dim)
        [Eigva, eigve] = H0.eigenstates()
        omega_a = Eigva[1] - Eigva[0]
        print((deltalambda - omega_a)**2., deltalambda, omega_a)

        return ((deltalambda - omega_a)/deltalambda)**2.


    def optimize_Efield(self, Deltalambda, eps_guess=0.):
        """
        Optimization routine to find the electric field based on the measured Delta omega and a guess for the field
        """
        from scipy.optimize import minimize

        res = minimize(self.find_Efield, x0=eps_guess, args=(self.molecule['B'], self.molecule['D'], self.dim, Deltalambda))
        print(res.x)

        return res.x

    @staticmethod
    def optimize_EfieldSpectrum(x, B, D, P, tau, reference_spectrum):
        """
        Obtain the efield from optimization of least square difference between spectra
        """
        eps = x[0]
        # HERE IS AN ERROR WITH RESPECT TO SELF
        c_sq = opt_2levelImpactEfield_interferometry(B, eps, D, P, tau, reference_spectrum)
        #c_sq = run_2levelImpactEfield_interferometry(eps, D, P, tau, reference_spectrum)
        self.inter = c_sq
        self.get_spectrum()

        lsq = ls(self.spectrum, reference_spectrum)

        return lsq


    def optimize_Efield2levelImpact(self, reference_spectrum, tau, P, eps0=0.):
        """
        Uses a least square fit to experimental or simulated spectrum to obtain the electric field strenth in the 2level impact approximation

        Args:
        -----
            reference spectrum: array of floats
                The 'correct' spectrum that we want to relate the field strength to.

            tau: array of floats
                The array of time delays used for the interferometry run

            P: float
                The pulse fluence used for the interferometry run

            eps0: float
                A first guess for the electric field. Default to zero.
        """
        from Utility import least_square_fitting as ls
        from scipy.optimize import minimize as mini

        #lsq = mini(optimize_EfieldSpectrum, x0=eps0, args=(self, self.D, P, tau, reference_spectrum))
        res = mini(self.optimize_EfieldSpectrum, x0=eps0, args=(self.molecule['B'], self.molecule['D'], P, tau, reference_spectrum))

        #c_sq = self.run_2levelImpactEfield_interferometry(eps0, self.D, P, tau)
        #self.inter = c_sq
        #self.get_spectrum()

        #lsq = ls(self.spectrum, reference_spectrum)

        return res.x


    def show_H0(self, scaleB=False):
        """
        """
        from Utility import H0
        
        if scaleB:
            B = 1.
        else:
            B = self.molecule.B
        H0 = H0(B, self.dim, self.molecule.a)
        return H0



class ImpactInterferometry(Interferometry):
    """
    Interferometry in the impact approximation
    """

    #def __init__(self, Ps=None, ts=None, B=1., istate=None,  mstate=0, dim=3, a=0., Name=None, odd=False, mrep=False) -> None:
    def __init__(self, Molecule:dict|str, Pulsepara:dict, istate=None,  mstate=0, dim=3, Name=None, odd=False, mrep=False) -> None:
        """
        Text
        """
        from Utility import ImpactPulses as IP
        #super().__init__(istate, mstate, B, dim, a, Name, odd, mrep)
        try:
            super().__init__(Molecule, istate, mstate, dim, Name, odd, mrep)
            self.Pulses = IP(Pulsepara)
        except Exception as e:
            print("An exception vas raised: {}".format(e))

    def initialize_pulses(self, Pulsepara) -> None:
        """
        Sets the pulses for the interferometry instance
            
            Args: 
            -----
                Pulsepara: dictionary
                    Contains Ps, taus
                    Ps: Array of floats 
                        Contains the pulse strengths 
                    taus: Array of floats
                        The time delays
        """
        #import Utility as Ut
        from Utility import ImpactPulses as IP

        self.Pulses = IP(Pulsepara)

        #Pulses = Ut.Pulses(Pulsepara['Ps'], Pulsepara['taus'])
        #self.Pulses = Pulses


    def set_pulses(self, Ps, taus) -> None:
        """
        Sets the pulses for the interferometry instance
            
            Args: 
            -----
                Ps: Array of floats 
                    Contains the pulse strengths 
                tauss: Array of floats
                    The time delays
        """
        #import Utility as Ut
        from Utility import ImpactPulses as IP
        
        Pulsepara = {'Ps': Ps, 'taus': taus}
        self.Pulses = IP(Pulsepara)

        #Pulses = Ut.Pulses(Ps, taus)
        #self.Pulses = Pulses

    def set_EvolutionOperator(self, P1, P2, tau1, tau2) -> None:
        """
        Sets the Evolution operator based on the pulse strengths and time delay
        """
        import Utility as Ut
        U1 = Ut.FullEvolutionOperator(P1, self.molecule.B, tau1, dim=self.dim, a=self.molecule.a, odd=self.odd, mrep=self.mrep)
        U2 = Ut.FullEvolutionOperator(P2, self.molecule.B, tau2, dim=self.dim, a=self.molecule.a, odd=self.odd, mrep=self.mrep)
        self.U = U2.U * U1.U

    def run_interferometry(self):
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        from qutip import ket2dm, basis
        from Utility import Proj
        lt = len(self.Pulses.taus)
        inter = np.zeros(lt)
        Op = ket2dm(basis(self.dim, self.mstate))
        for i in range(lt):
            print("Running interferometry, step {} of {}".format(i+1, lt))
            #self.set_pulses(self.Pulses.Ps, [0., self.Pulses.taus[i]])
            self.set_EvolutionOperator(self.Pulses.Ps[0], self.Pulses.Ps[1], 0., self.Pulses.taus[i])
            fstate = self.U * self.istate
            dm = ket2dm(fstate)
            inter[i] = Proj(Op,dm) # Check real not complex!!!
            #fspop  = fstate.extract_states(self.mstate)
            #inter[i] = abs(fspop.overlap(fspop))**2. 
            
        
        self.inter = inter

    def print_info(self):
        """
        """
        super().print_info()
        print()
        self.Pulses.print_pulse_info()

class ImpactEfieldInterferometry(ImpactInterferometry):
    """
    Interferometry in the impact approximation
    """


    #def __init__(self, Ps=None, ts=None, B=1., eps=0., D=0., istate=None,  mstate=0, dim=3, a=0., Name=None, odd=False, mrep=False, alpha=None, beta=None, pol=None) -> None:
    def __init__(self, Molecule, Pulsepara, Fieldpara, istate=None,  mstate=0, dim=3, Name=None, odd=False, mrep=False) -> None:
        """
        Text
        """
        from Utility import ImpactPulses as IP, StaticEfield as EF
        try:
            # Setting molecule and the impact interferometry pulses, initial state, dimension etc 
            super().__init__(Molecule, Pulsepara, istate, mstate, dim, Name, odd, mrep)
            if hasattr(self.Pulses, 'alpha'):
                print("alpha", self.Pulses.alpha)
            else:
                print("No alpha")

            #if "alpha" not in Pulsepara:
            #    self.Pulses.alpha = None
            #else:
            #    self.Pulses.alpha = Pulsepara['alpha']
            #    #self.Pulses.alpha = np.zeros(len(self.Pulses.Ps))
            #if "beta" not in Pulsepara:
            #    self.Pulses.beta = None
            #else:
            #    self.Pulses.beta = Pulsepara['beta']
                #self.Pulses.beta = np.zeros(len(self.Pulses.Ps))
            # Setting the static electric field
            self.EF = EF(Fieldpara)
            #self.set_EvolutionOperator()

        except Exception as e:
            print("An exception occured in initialization: {}".format(e))
        #if alpha is None:
        #    super().__init__(istate, mstate, B, dim, a, Name, odd, mrep)
        #else:
        #    super().__init__(istate, mstate, B, dim, a, Name, odd, mrep=True)

        #self.Ps       = Ps
        #self.ts       = ts
        #self.eps      = eps
        #self.alpha    = alpha
        #self.beta     = beta
        #self.D        = D
        #self.tau      = ts
        #self.name     = Name
        #self.pol      = pol


        #if (Ps is not None) and (ts is not None):
        #    self.set_pulses(Ps, ts, self.alpha, self.beta)

                
    
    def initialize_Efield(self, para:dict):
        """
        Sets the static electric field we want to measure.
        """
        from Utility import Static_Efield as EF
        self.Efield = EF(para)

    def set_pulses(self, Ps, ts, alpha=None, beta=None) -> None:
        """
        Sets the pulses for the interferometry instance
            
            Args: 
                Ps: Array of floats 
                    Contains the pulse strengths 
                ts: Array of floats
                    The time delays
                alpha: Array of floats
                    The polarization angles of the pulses w.r.t. the z-axis in the lab frame
                beta: Array of floats
                    The polarization angles of the pulses w.r.t. the x-y-axis in the lab frame
        """
        from Utility import ImpactPulses as IP
        if self.Pulses.alpha is None:
            Pulsepara = {'Ps': self.Pulses.Ps, 'taus': self.Pulses.taus}
        else:
            Pulsepara = {'Ps': self.Pulses.Ps, 'taus': self.Pulses.taus, 'alpha': self.Pulses.alpha, 'beta': self.Pulses.beta}

        #if self.Pulses.alpha.any() == None:
        #    alpha = np.zeros(self.Pulses.nP)
        #if self.Pulses.beta.any() == None:
        #    beta = np.zeros(self.Pulses.nP)
        #angles = {'alpha': alpha, 'beta': beta}
        #Pulses = Ut.Pulses(Ps, ts, angles)
        self.Pulses = IP(Pulsepara)
        #if pol is not None:
        #    self.Pulses.pol = pol

    def set_alpha(self, alpha, index) -> None:
        if (index == 0) or index == 1:
            self.alpha[index] = alpha
        else:
            raise ValueError("index must be either 0 or 1")
    

    def set_beta(self, beta, index) -> None:
        if (index == 0) or index == 1:
            self.beta[index] = beta
        else:
            raise ValueError("index must be either 0 or 1")

    
    def set_EvolutionOperator(self, time_ind) -> None:
        """
        Sets the Evolution operator based on the pulse strengths and time delay
        """
        from Utility import FullEfieldEvolutionOperator,  U2U1
        from qutip import Qobj, qeye
        #if hasattr(self.Pulses, 'pol'):
        #if self.mrep:
        #    n = (self.dim + 1)**2.
        #else:
        #    n = self.dim
        #Uhold = Qobj(qeye(n))
        print(time_ind, self.Pulses.taus[time_ind])
        if hasattr(self.Pulses.alpha, 'array') and hasattr(self.Pulses.beta, 'array'):
            Uhold = Qobj(qeye((self.dim+1)**2))
            FEE1 =FullEfieldEvolutionOperator(self.Pulses.Ps[0], self.molecule.B, 0., self.EF.eps0, self.molecule.D, (self.dim+1)**2, self.molecule.a, self.Name, self.odd, self.mrep, True, self.Pulses.alpha[0], self.Pulses.beta[0])  
            FEE2 =FullEfieldEvolutionOperator(self.Pulses.Ps[1], self.molecule.B, self.Pulses.taus[ime_ind], self.EF.eps0, self.molecule.D, (self.dim+1)**2, self.molecule.a, self.Name, self.odd, self.mrep, True, self.Pulses.alpha[1], self.Pulses.beta[1])  

            #for i in range(self.Pulses.nP):
                #FEE =FullEfieldEvolutionOperator(self.Pulses.Ps[i], self.molecule.B, self.Pulses.taus[i], self.EF.eps0, self.molecule.D, (self.dim+1)**2, self.molecule.a, self.Name, self.odd, self.mrep, True, self.Pulses.alpha[i], self.Pulses.beta[i])  
            Uhold = U2U1(FEE2.U, FEE1.U) 
            self.U = Uhold
            #U1 = FullEfieldEvolutionOperator(self.Pulses.P[0], self.B, self.Pulses.t[0], self.eps, self.D, self.dim, self.a, self.name, self.odd, self.mrep, True, self.alpha[0], self.beta[0])
            #U2 = FullEfieldEvolutionOperator(self.Pulses.P[1], self.B, self.Pulses.t[1], self.eps, self.D, self.dim, self.a, self.name, self.odd, self.mrep, True, self.alpha[1], self.beta[1])
        else:
            Uhold = Qobj(qeye(self.dim))
            FEE1 = FullEfieldEvolutionOperator(self.Pulses.Ps[0], self.molecule.B, 0., self.EF.eps0, self.molecule.D, dim=self.dim, a=self.molecule.a, name=self.Name, odd=self.odd, mrep=self.mrep) 

            FEE2 = FullEfieldEvolutionOperator(self.Pulses.Ps[1], self.molecule.B, self.Pulses.taus[time_ind], self.EF.eps0, self.molecule.D, dim=self.dim, a=self.molecule.a, name=self.Name, odd=self.odd, mrep=self.mrep) 
            #for i in range(self.Pulses.nP):
            #    FEE = FullEfieldEvolutionOperator(self.Pulses.Ps[i], self.molecule.B, self.Pulses.taus[i], self.EF.eps0, self.molecule.D, dim=self.dim, a=self.molecule.a, name=self.Name, odd=self.odd, mrep=self.mrep) 
            Uhold = U2U1(FEE2.U, FEE1.U)
            self.U = Uhold
 
            #U1 = FullEfieldEvolutionOperator(self.Pulses.P[0], self.B, self.Pulses.t[0], self.eps, self.D, dim=self.dim, a=self.a, name=self.name, odd=self.odd, mrep=self.mrep)
            #U2 = FullEfieldEvolutionOperator(self.Pulses.P[1], self.B, self.Pulses.t[1], self.eps, self.D, dim=self.dim, a=self.a, name=self.name, odd=self.odd, mrep=self.mrep)
        #self.U = U2U1(U2.U, U1.U) #U2.U * U1.U


    def run_interferometry(self):
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        from qutip import ket2dm, basis, qeye
        from Utility import Proj, index
        lt = len(self.Pulses.taus)
        inter = np.zeros(lt)
        if self.mrep:
            Op = qeye((self.dim+1)**2)
            for m in range(-self.mstate, self.mstate+1):
                Op += ket2dm(basis((self.dim + 1)**2, index(self.mstate, m)))
            Op -= qeye((self.dim+1)**2)
        else:
            Op = ket2dm(basis(self.dim, self.mstate))
        for i in range(lt):
            #self.set_pulses(self.Pulses.Ps, [0., self.Pulses.taus[i]], self.Pulses.pol, self.Pulses.alpha)
            self.set_EvolutionOperator(i)
            fstate = self.U * self.istate
            dm = ket2dm(fstate)
            inter[i] = Proj(Op,dm)
            print(i, inter[i])
            #fspop  = fstate.extract_states(self.mstate)
            #inter[i] = abs(fspop.overlap(fspop))**2. 
            
        
        self.inter = inter


    def find_efield(x, B, D, dim, omega, mrep=False):
        """
        Find the electric field from the molecular parameters and the measured value of omega
        """

        eps = x[0]

        if mrep:
            H = Ut.H0_m(B, dim, full=True) + Ut.Hdip_m(eps, D, dim)
            [Eigva, Eigv] = H.eigenstates()
            omega_a = Eigva[3] - Eigva[0]
        else:
            H = Ut.H0(B, dim, full=True) + Ut.Hdip(eps, D, dim)
            [Eigva, Eigv] = H.eigenstates()
            omega_a = Eigva[1] - Eigva[0]

        return ((omega - omega_a) / omega)**2.


    def optimize_efield(eps0=1.e-6, omega=0.):
        """
        Find the efield by optimization from the initial guess for the efield and the measured energy difference
        """
        from scipy.optimize import minimize as mini

        res = mini(find_efield, x0=eps0, args=(self.B, self.D, self.dim, omega*self.B))
        eps = res.x[0]

        return eps


    def find_eps0(omin=None, omax=None, p0=None, omi=200.):
        """
        """
        omega_m = self.find_omega0(omin, omax, p0, omi) * self.B
        #Using the two level approximation to obtain an estimate for eps0

        eps0 = math.sqrt(3.* (0.25*(omega_m)**2. - self.B**2.)) / self.D
        
        return eps0



class ImpactTensorInterferometry(RotationTensor, ImpactInterferometry):
    """
    Interferometry in the impact approximation with coupled molecular ions
    """

    def __init__(self, Molpara, Pulsepara, Molpara2=None, dim=3, istate=None, mstate=0, odd=False, mrep=False, name=None, set_Hams=False):
        """
        """
        from qutip import basis, tensor
        print("ImpactTensorInterferomerty UNDER DEVELOPMENT!!!")
        try:
            RotationTensor.__init__(self, Molpara, Molpara2, dim=dim, istate1=istate, name=name)
            self.initialize_pulses(Pulsepara)
            if istate is None:
                self.istate = tensor(basis(dim, 0), basis(dim, 0))
            else:
                self.istate = istate
            if mstate is None:
                self.mstate = 0
            else:
                self.mstate = mstate
            if set_Hams:
                self.set_Ham()
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError {}:".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))
            

    def print_istate(self, state):
        """
        """
        print("Initial state (psi):")
        print("--------------")
        print(self.istate)
        print()
        print("Initial state (phi):")
        print("--------------")
        print(self.CX_psi(self.istate))
        print()


    def print_info(self):
        """
        """
        super().print_info(self)
        self.Pulses.print_info()
        self.print_istate(self.istate)


    def run_interferometry(self):
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        from qutip import ket2dm, basis, ptrace
        from Utility import Proj
        lt        = len(self.Pulses.taus)
        inter_psi = np.zeros(lt)
        inter_phi = np.zeros(lt)
        #Op        = ket2dm(tensor(basis(self.dim[0], self.mstate), basis(self.dim[0], self.mstate)))
        Op1        = ket2dm(tensor(basis(self.dim[0], 0), basis(self.dim[0], 0)))
        Op2        = ket2dm(tensor(basis(self.dim[0], 0), basis(self.dim[0], 1)))
        Op = Op1 + Op2
        self.set_impact_U(self.Pulses.taus[50], self.Pulses.Ps[0], self.Pulses.Ps[1])
        U = self.Ui12 * self.Uf * self.Ui11

        # First run, psi
        print("NOT IN WORKING CONDITION!!!")
        for i in range(lt):
            self.set_impact_U(self.Pulses.taus[i], self.Pulses.Ps[0], self.Pulses.Ps[1])
            #self.set_impact_U(self.Pulses.taus[i], 0., 0.)
            #self.set_pulses(self.Pulses.Ps, [0., self.Pulses.taus[i]])
            #self.set_EvolutionOperator()
            #fstate_psi   = self.U * self.istate
            U = self.Ui12 * self.Uf * self.Ui11
            fstate_psi   = U * self.istate
            #fstate_phi   = self.U * self.CX_psi(istate)
            fstate_phi   = U * self.CX_psi(self.istate)
            dm_psi       = ket2dm(fstate_psi)
            dm_phi       = ket2dm(fstate_phi)
            inter_psi[i] = Proj(Op,dm_psi)
            inter_phi[i] = Proj(Op,dm_phi)
            
        
        self.inter_psi = inter_psi
        self.inter_phi = inter_phi




class PulsesInterferometry(Interferometry):
    """
    Interferometry taking the full effect o the pulses into account
    """

    def __init__(self, Molpara, Pulsepara=None, istate=None, mstate=0, dim=3, Name=None, odd=False, mrep=False, custom_para=False) -> None:
        """
        Initialization routine for a two orr-resonance fs-pulse interferometry run

            Args:
            -----

            Pulsepara: dict
                Contains the pulse paramters: I0, sigma, t0

            time: array of floats
                Time grid for the interferometer

            molpara: dict
                Contains the molecular parameters: B, a and Da

            istate: Qobj
                The initial state

            mstate: int
                The state on which to measure the final population

            dim: int
                The dimenstion of the basis

            Name: str
                Name of the instance

            odd: Boolean
                Whether to use only the odd states

            mrep: Boolean
                Whether to include all mstates (True) o rnot (False)
        """
        from Utility import GaussPulses as GP
        try:
            super().__init__(Molpara, istate, mstate, dim, Name, odd, mrep, custom_para=None) 
            if Pulsepara is not None:
                self.Pulses = GP(Pulsepara)
            else:
                self.Pulses = None
        except KeyError as e:
            raise Exception(e)
        
    def set_Hams(self):
        """
        """
        from Utility import H0, H1_I0, H0_m

        if self.mrep:
            print("WARNING: Only full representation implemented for m representation yet!")
            H0 = H0_m(self.molecule.B, self.dim, self.molecule.a, self.odd, full=True)
            HI = H1_I0(self.molecule.Da, self.dim, self.odd, full=True)
        else:
            H0 = H0(self.molecule.B, self.dim, self.molecule.a)
            HI = H1_I0(self.molecule.Da, self.dim)

        return H0, HI
        
    def run_interferometry(self, options=None, split=True) -> None:
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        from qutip import ket2dm, basis, Options
        from Utility import H0, H1_I0, Proj, FreeEvolutionOperator as FE
        
        # Set the free and interaction terms of the Hamiltonian
        H0, HI = self.set_Hams()
        
        # Define the measurment operator and the interferometry array
        ltau = len(self.Pulses.taus)
        Op = ket2dm(basis(self.dim, self.mstate))
        inter = np.zeros(ltau)

        # Options to the numerical qutip solver
        if options is None:
            options = Options(store_final_state=True)
        else:
            if options.store_final_state == False:
                print("Warning, must store final state in order to record final state population. Setting to True!")
                options.store_final_state = True

        # Run the interferometry scheme. We may chose to split the time grid for eficeiency
        if split:
            exp_fact = 5.
            for i in range(ltau):
                print("Running interferometry, step {} of {}".format(i+1, ltau))
                if self.Pulses.taus[i] > 2. * exp_fact * self.Pulses.sigma:
                    free_time = self.Pulses.taus[i] - 2. * exp_fact * self.Pulses.sigma
                    FreeU = FE(self.molecule.B, free_time, self.dim, self.molecule.a)
                    #time1 = np.arange(self.Pulses.tmin, self.Pulses.tmax, self.Pulses.dt)
                    time1 = np.arange(self.Pulses.tmin, exp_fact * self.Pulses.sigma, self.Pulses.dt)
                    time2 = np.arange(exp_fact * self.Pulses.sigma + free_time, self.Pulses.taus[i] - self.Pulses.tmin, self.Pulses.dt)
                    #output = mesolve([H0, [HI, double_Gauss_me]], self.istate, time, e_ops=[Op], options=options, args={'t0': 0.,'I01': self.Pulses.I0s[0], 'I02': self.Pulses.I0s[1],'sigma': self.Pulses.sigma, 'tau': self.Pulses.taus[i]})
                    #output1 = sesolve([H0, [HI, Gauss_me]], self.istate, time1, e_ops=[Op], options=options, args={'t0': 0.,'I0': self.Pulses.I0s[0], 'sigma': self.Pulses.sigma})
                    #fstate1 = output1.final_state
                    fstate1 = self.run_Gauss(H0, HI, self.istate, time1, 0., 0, Op, options)
                    freestate = FreeU.Uf * fstate1
                    #output = sesolve([H0, [HI, Gauss_me]], freestate, time2, e_ops=[Op], options=options, args={'t0': self.Pulses.taus[i], 'I0': self.Pulses.I0s[1], 'sigma': self.Pulses.sigma})
                    fstate2 = self.run_Gauss(H0, HI, freestate, time2, self.Pulses.taus[i], 1, Op, options)
                    #dm = ket2dm(output.final_state)
                    dm = ket2dm(fstate2)
                    inter[i] = Proj(Op, dm)
                else:
                    #time = np.arange(self.Pulses.tmin, self.Pulses.tmax, self.Pulses.dt)
                    time = np.arange(self.Pulses.tmin, self.Pulses.taus[i] - self.Pulses.tmin, self.Pulses.dt)
                    #output = sesolve([H0, [HI, double_Gauss_me]], self.istate, time, e_ops=[Op], options=options, args={'t0': 0.,'I01': self.Pulses.I0s[0], 'I02': self.Pulses.I0s[1],'sigma': self.Pulses.sigma, 'tau': self.Pulses.taus[i]})

                    fstate = self.run_doubleGauss(H0, HI, time, i, Op, options)
                    #dm = ket2dm(output.final_state)
                    dm = ket2dm(fstate)
                    inter[i] = Proj(Op, dm)
        else:
            time = np.arange(self.Pulses.tmin, self.Pulses.tmax, self.Pulses.dt)
            for i in range(ltau):
                print("Running interferometry, step {} of {}".format(i+1, ltau))
                #ime = np.arange(self.Pulses.tmin, self.Pulses.tmax, self.Pulses.dt)
                #utput = sesolve([H0, [HI, double_Gauss_me]], self.istate, time, e_ops=[Op], options=options, args={'t0': 0.,'I01': self.Pulses.I0s[0], 'I02': self.Pulses.I0s[1],'sigma': self.Pulses.sigma, 'tau': self.Pulses.taus[i]})
                fstate = self.run_doubleGauss(H0, HI, time, i, Op, options)
                #dm = ket2dm(output.final_state)
                dm = ket2dm(fstate)
                inter[i] = Proj(Op, dm)
        self.inter = inter

    def run_doubleGauss(self, H0, HI, time, time_index, Op, options):
        """
        """
        from Utility import double_Gauss_me

        output = sesolve([H0, [HI, double_Gauss_me]], self.istate, time, e_ops=[Op], options=options, args={'t0': 0.,'I01': self.Pulses.I0s[0], 'I02': self.Pulses.I0s[1],'sigma': self.Pulses.sigma, 'tau': self.Pulses.taus[time_index]})
        
        return output.final_state


    def run_Gauss(self, H0, HI, istate, time, t0, pulse_index, Op, options):
        """
        """
        from Utility import Gauss_me 
        
        output = sesolve([H0, [HI, Gauss_me]], istate, time, e_ops=[Op], options=options, args={'t0': t0,'I0': self.Pulses.I0s[pulse_index], 'sigma': self.Pulses.sigma})

        return output.final_state



    def print_info(self):
        """
        """
        super().print_info()
        self.Pulses.print_info()
        #self.print_istate(self.istate)



class PulsesEfieldInterferometry(PulsesInterferometry):
    """
    Interferometry taking the full effect o the pulses into account
    """

    def __init__(self, Molecule, Pulsepara=None,  Fieldpara=None, istate=None, mstate=0, dim=3, Name=None, mrep=False, custom_para=False) -> None:
        """
        Initialization routine for a two orr-resonance fs-pulse interferometry run

            Args:
            -----

            Pulsepara: dict
                Contains the pulse paramters, I0, sigma, t0

            time: array of floats
                Time grid for the interferometer

            molpara: dict
                Contains the molecular parameters, B, dip and Da

            istate: Qobj
                The initial state

            mstate: int
                The state on which to measure the final population

            dim: int
                The dimenstion of the basis

            Name: str
                Name of the instance
        """
        from Utility import StaticEfield as SE
        try:
            super().__init__(Molecule, Pulsepara, istate, mstate, dim, Name, mrep, custom_para)
            self.SE = SE(Fieldpara)
        except KeyError as e:
            raise Exception(e)
   

    def set_Hams(self):
        """
        """
        from Utility import H0, H0_m, H1_I0, H1_I0_m, H_dip, H_dipm

        if self.mrep:
            H0 = H0_m(self.molecule.B, self.dim, self.molecule.a, full=True) + H_dipm(self.SE.eps0, self.molecule.D, self.dim)
            HI = H1_I0_m(self.molecule.Da, self.dim, full=True) 
        else:
            H0 = H0(self.molecule.B, self.dim, self.molecule.a, full=True) + H_dip(self.SE.eps0, self.molecule.D, self.dim)
            HI = H1_I0(self.molecule.Da, self.dim, full=True)

        return H0, HI
     

    def run_interferometry(self, options=None) -> None:
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        #from qutip import mesolve, ket2dm, basis, Options
        from qutip import ket2dm, basis, Options
        #from Utility import H0, H1_I0, H_dip, double_Gauss_me, Proj
        from Utility import H0, H1_I0, H_dip, Proj
        
        #if self.mrep:
        #    H0 = H0_m(self.molecule.B, self.dim, self.molecule.a, full=True) + H_dipm(self.SE.eps0, self.molecule.D, self.dim)
        #    HI = H1_I0_m(self.molecule.Da, self.dim, full=True) 
        #else:
        #    H0 = H0(self.molecule.B, self.dim, self.molecule.a, full=True) + H_dip(self.SE.eps0, self.molecule.D, self.dim)
        #    HI = H1_I0(self.molecule.Da, self.dim, full=True)

        # Set the free and interaction terms of the Hamiltonian
        H0, HI = self.set_Hams()
        
        # Define the measurment operator and the interferometry array
        ltau = len(self.Pulses.taus)
        Op = ket2dm(basis(self.dim, self.mstate))
        inter = np.zeros(ltau)

        # Options to the numerical qutip solver
        if options is None:
            options = Options(store_final_state=True)
        else:
            if options.store_final_state == False:
                print("Warning, must store final state in order to record final state population. Setting to True!")
                options.store_final_state = True

        # Run the interferometry scheme. We may chose to split the time grid for eficeiency
        time = np.arange(self.Pulses.tmin, self.Pulses.tmax, self.Pulses.dt)
        for i in range(ltau):
            print("Running interferometry, step {} of {}".format(i+1, ltau))
            fstate = self.run_doubleGauss(H0, HI, time, i, Op, options)
            #output = mesolve([H0, [HI, double_Gauss_me]], self.istate, time, e_ops=[Op], options=options, args={'t0': 0.,'I01': self.Pulses.I0s[0], 'I02': self.Pulses.I0s[1],'sigma': self.Pulses.sigma, 'tau': self.Pulses.taus[i]})

            #dm = ket2dm(output.final_state)
            dm = ket2dm(fstate)
            inter[i] = Proj(Op, dm)
        self.inter = inter



