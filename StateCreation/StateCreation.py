"""Module containing state creation classes and methods in the impulse approximation and without the approximation"""

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


"""
Some useful methods
"""

def set_internal_operators(Ps, taus, B, dim):
    from Utility import EvolutionOperators as EOs
    from Utility import ImpactPulses as IP
    Pulses = IP(Ps, taus)
    return EOs(Pulses, B, dim)

"""
Classes
"""

class StateCreation():
    """
    The main state creation class
    """

    def __init__(self, Molecule:str, Pulsepara:dict, istate, goal_state, Name:str=None, approximation:str=None, measure="overlap") -> None:
        """
        """
        self.Molecule   = pm.set_MolParams(Molecule)
        self.Pulsepara  = Pulsepara
        self.istate     = istate
        self.goal_state = goal_state
        self.dim        = istate.dims[0][0]
        self.measure    = measure
        if Name == None:
            self.Name = ""
        else:
            self.Name = Name
        if approximation == None:
            self.appr = ""
        else:
            self.appr = approximation

        if approximation == 'impact':
            try:
                if len(Pulsepara['Ps']) == len(Pulsepara['taus']):
                    self.Ps  = Pulsepara['Ps']
                    self.taus = Pulsepara['taus']
                else:
                    print("Number of pulses and time delays must be the same. Got {} pulses and {} time delays".format(len(Pulsepara['Ps']), len(Pulsepara['taus'])))
            except KeyError as e:
                raise Exception(e)
        else:
            print("Implement full pulse initialization!!!")


    def set_operators(self):
        from Utility import EvolutionOperators as EOs
        from Utility import ImpactPulses as IP
        Pulses = IP(self.Ps, self.taus)
        self.EOS = EOs(Pulses, self.Molecule['B'], self.dim)

    def update_operators(self, Ps=None, taus=None):
        from Utility import EvolutionOperators as EOs
        from Utility import ImpactPulses as IP
        if Ps is None:
            Ps_new = self.Ps
        else:
            Ps_new = Ps
        if taus is None:
            taus_new = self.taus
        else:
            taus_new = taus

        Pulses = IP(Ps_new, taus_new)
        self.EOS = EOs(Pulses, self.Molecule['B'], self.dim)




    @staticmethod
    def opt_P_tau_impact_overlap(x, B, dim, istate, goalstate):
        """
        Minimizer for state overlap in the impact approximation
        """
        # Seprate out the pulse strengths and time delays
        mid = len(x)//2
        Ps = x[:mid]
        taus = x[mid:]

        # Make sure that both are positive
        l2 = 0.
        for P in Ps:
            if P < 0.:
                l2 = 10000.
        for tau in taus:
            if tau < 0.:
                l2 = 10000.

        Op = set_internal_operators(Ps, taus, B, dim) # Set the impact and free operators 
        fstate = Op.U * istate                        # Propagate the initial state
        ov = goalstate.overlap(fstate)                # Get the overlap
        l1 = abs(1. - abs(ov))**2.                    # Error function
        return l1 + l2 

    def run_sequence(self):

        fstate = self.EOS.U * self.istate
        return fstate

    def optimize_pulse_sequence(self):
        """
        """
        from scipy.optimize import minimize as mini
        
        if self.appr == "impact":
            self.set_operators()
            res = mini(self.opt_P_tau_impact_overlap, x0=np.concatenate((self.Ps, self.taus)), \
                    args=(self.Molecule['B'], self.dim, self.istate, self.goal_state))
            return res
            #self.run_optimize_impact()
        else:
            print("Unknown approximation {]".format(self.appr))


class StateCreation2D(StateCreation):
    """
    State cration in a 2-level basis (|j=0,m=0> and |j=2,m=0>)
    """

    def __init__(self, Molecule:dict, Pulsepara:dict, istate, goal_state, approximation:str) -> None:
        """
        """
        # Check that istate and goal_state are Qobj, either density matrices or kets and 2-level
        if ((istate.type == 'ket') and (goal_state.type == 'ket') or (istate.type == 'oper') and (goal_state.type == 'oper')):
            if ((istate.type == 'ket') and (istate.dims[0][0] == 2) or (istate.type == 'oper') and (istate.dims[0][0] == 2) and (istate.dims[1][1] == 2)):
                super().__init__(Molecule, Pulsepara, istate, goal_state, "2-level", approximation)
        else:
            raise Exception





