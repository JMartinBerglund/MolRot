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

class Interferometry():
    """
    The main interferometry class
    """

    def __init__(self, istate, mstate=0, B=1., dim=3, Name=None) -> None:
        """
        Initializes an instance of the interferometer
            
            Args:
                istate: Qobj (ket or density matrix)
                    Initial state for the interferometer. 
                
                mstate : int
                    Which final state to measure the interferogram on
                
                B: float
                    The rotational constant

                dim: int
                    The dimansion of the basis used

                Name: str
                    Name of the interferomoeter intance

        """
        self.B        = B
        self.istate   = istate
        self.mstate   = mstate
        self.dim      = dim
        self.Name     = Name
        self.U        = None
        self.inter    = None
        self.tau      = None
        self.spectrum = None
        self.omega    = None

    def run_interferometry(self) -> None:
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        pass

    def get_spectrum(self, shift=False) -> None:
        """
        Obtains the Fourier spectrum of the interferogram
        """
        from scipy.fft import fft, fftfreq, fftshift
        # Check that the interferogram and delay arrays exist
        # and take the Fourier transform if they do
        if (self.tau is not None) and (self.inter is not None):
            ntau = len(self.tau)
            dtau = self.tau[1] - self.tau[0]
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
        
    def plot_interferogram(self, title=None) -> None:
        import matplotlib.pyplot as plt
        if (self.tau is not None) and (self.inter is not None):
            # Prepare plot
            plt.plot(self.tau, self.inter)

            # Set axis labels
            plt.xlabel('$\\tau$')
            plt.ylabel('Interferogram')

            # Set the title
            if title is not None:
                plt.title(title)

            # Plot the graph
            plt.plot()
        else:
            print("Tried to plot the Fourier spectrum. Can't do that if either the spectrum or omega grid don' t exist.")
            raise Exception()

    def plot_spectrum(self, title=None, scaleTrot=False, setTicks=False) -> None:
        import matplotlib.pyplot as plt
        if (self.omega is not None) and (self.spectrum is not None):
            # Prepare plot
            if scaleTrot is True:
                Trot = 2.*math.pi*(self.B)**(-1.)
                plt.plot(self.omega * Trot, self.spectrum)
            else:
                plt.plot(self.omega, self.spectrum)


            # Set axis labels 
            if scaleTrot:
                plt.xlabel('$\omega$ ($T_{rot}^{-1}$)' )
                if setTicks:
                    plt.xticks(np.arange(self.omega[0]*Trot, self.omega[-1]*Trot, step=2))
            else:
                plt.xlabel('$\omega$ (au)')

                    
            plt.ylabel('Spectrum')

            # Set the title
            if title is not None:
                plt.title(title)

            # Plot the graph
            plt.plot()
            plt.show()
        else:
            print("Tried to plot the Fourier spectrum. Can't do that if either the spectrum or omega grid don' t exist.")
            raise Exception()


class ImpactInterferometry(Interferometry):
    """
    Interferometry in the impact approximation
    """

    def __init__(self, Ps=None, ts=None, B=1., istate=None,  mstate=0, dim=3, Name=None) -> None:
        """
        Text
        """
        super().__init__(istate, mstate, B, dim, Name)
        self.Ps       = Ps
        self.tau      = ts
        


    def set_pulses(self, Ps, ts) -> None:
        """
        Sets the pulses for the interferometry instance
            
            Args: 
                Ps: Array of floats 
                    Contains the pulse strengths 
                ts: Array of floats
                    The time delays
        """
        import Utility as Ut
        Pulses = Ut.Pulses(Ps, ts)
        self.Pulses = Pulses

    def set_EvolutionOperator(self) -> None:
        """
        Sets the Evolution operator based on the pulse strengths and time delay
        """
        import Utility as Ut
        U1 = Ut.FullEvolutionOperator(self.Pulses.P[0], self.B, self.Pulses.t[0], dim=self.dim)
        U2 = Ut.FullEvolutionOperator(self.Pulses.P[1], self.B, self.Pulses.t[1], dim=self.dim)
        self.U = U2.U * U1.U

    def run_interferometry(self):
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        from qutip import ket2dm, basis
        from Utility import Proj
        lt = len(self.tau)
        inter = np.zeros(lt)
        Op = ket2dm(basis(self.dim, self.mstate))
        for i in range(lt):
            self.set_pulses(self.Ps, [0., self.tau[i]])
            self.set_EvolutionOperator()
            fstate = self.U * self.istate
            dm = ket2dm(fstate)
            inter[i] = Proj(Op,dm)
            #fspop  = fstate.extract_states(self.mstate)
            #inter[i] = abs(fspop.overlap(fspop))**2. 
            
        
        self.inter = inter

class PulsesInterferometry(Interferometry):
    """
    Interferometry taking the full effect o the pulses into account
    """

    def __init__(self, Pulsepara=None, time=[], molpara=None, istate=None, mstate=0, dim=3, Name=None) -> None:
        """
        Initialization routine for a two orr-resonance fs-pulse interferometry run

            Args:
            -----

            Pulsepara: dict
                Contains the pulse paramters, I0, sigma, t0

            time: array of floats
                Time grid for the interferometer

            molpara: dict
                Contains the molecular parameters, B and Da

            istate: Qobj
                The initial state

            mstate: int
                The state on which to measure the final population

            dim: int
                The dimenstion of the basis

            Name: str
                Name of the instance
        """
        import Utility as Ut
        try:
            super().__init__(istate, mstate, molpara['B'], dim, Name)
        except KeyError as e:
            raise Exception(e)
        
        if Pulsepara is not None:
            try:
                self.I0    = [Pulsepara['I01'], Pulsepara['I02']]
                self.sigma = Pulsepara['sigma']
                self.t0    = Pulsepara['t0']
                self.tau   = Pulsepara['tau']
                self.time  = time
            except KeyError as e:
                raise Exception(e)
        else:
            self.I0    = None
            self.sigma = None
            self.t0    = None
            self.tau   = None
        if molpara is not None:
            try:
                self.Da = molpara['Da']
            except KeyError as e:
                raise Exception(e)
        else:
            self.Da = None

        #print(Ut.getP_int(self.Da, 2.*self.I0[0], self.t0, self.sigma))



    #def set_pulses(self, Pulses) -> None:
    #    """
    #    Sets the pulses for the interferometry instance
    #            Pulses: Array of pulses instance 
    #                Contains the pulse strengths and time delay
    #    """
    #    import Utility as Ut
    #    nP  = len(Pulses)
    #    self.Pulses  = []
    #    for i in range(nP):
    #        Pulse = Ut.Pulses(Pulses[i].P, Pulses[i].t)
    #        self.Pulses.append(Pulse)

    #def set_EvolutionOperators(self) -> None:
    #    """
    #    Text
    #    """
    #    import Utility as Ut
    #    nP  = len(self.Pulses)
    #    print("Number of pulses:", nP)
    #    self.Us  = []
    #    for i in range(nP):
    #        self.Us.append(Ut.EvolutionOperators(self.Pulses[i], self.B, dim=self.dim))

        
    def run_interferometry(self, options=None) -> None:
        """
        Runs an interferometry run and records final time delay dependent populatins
        """
        from qutip import mesolve, ket2dm, basis, Options
        from Utility import H0, H1_I0, double_Gauss_me, Proj
        H0 = H0(self.B, self.dim)
        HI = H1_I0(self.Da, self.dim)
        ltau = len(self.tau)
        Op = ket2dm(basis(self.dim, self.mstate))
        inter = np.zeros(ltau)
        if options is None:
            options = Options(store_final_state=True)
        else:
            if options.store_final_state == False:
                print("Warning, must store final state in order to record final state population. Setting to True!")
                options.store_final_state = True
        for i in range(ltau):
            output = mesolve([H0, [HI, double_Gauss_me]], self.istate, self.time, e_ops=[Op], options=options, args={'t0': self.t0,'I01': self.I0[0], 'I02': self.I0[1],'sigma': self.sigma, 'tau': self.tau[i]})

            dm = ket2dm(output.final_state)
            inter[i] = Proj(Op, dm)
        self.inter = inter
