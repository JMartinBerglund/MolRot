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

    def __init__(self, Molecule:str, istate, goal_state) -> None:
        """
        """
        self.Molecule = pm.sel_MolParams(Molecule)
        self.istate = istate
        self.goal_state = goal_state



