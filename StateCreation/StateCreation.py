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

class StateCreation():
    """
    The main state creation class
    """

    def __init__(self, Molecule:str, Pulsepara:dict, istate, goal_state, Name:str=None, approximation:str=None) -> None:
        """
        """
        self.Molecule   = pm.set_MolParams(Molecule)
        self.Pulsepara  = Pulsepara
        self.istate     = istate
        self.goal_state = goal_state
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


    def set_operators(self, Pulses:dict):
        pass


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





