"""Module containing molecular vibrarional classes and methods.
"""

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


class Vibration():
    """
    General class for handling molecular vibrations
    """

    def __init__(self, Molpara, dim=3, Vtype=None, Vname=None) -> None:
        """
        The generator for the molecular vibraions class.
        """
        try:
            self.check_molpara(Molpara, ["mu"])
            self.mu   = Molpara['mu']
            self.dim  = dim
            self.type = Vtype
            self.name = Vname
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError {}:".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))

    def check_molpara(self, para, keys=[]) -> None:
        """
        Checking the para file for inconsistencies.
        """
        for key in keys:
            if key not in para:
                raise KeyError("key {} not in para.".format(key))
            elif not isinstance(para[key], float):
                raise TypeError("{} needs to be of type float, got {}.".format(key, type(para[key])))
            elif para[key] < 0.:
                raise ValueError("{} needs to be a positive number, got {}".format(key, para[key]))

    def print_info(self):
        """
        Prints info about the instance to screen
        """
        if self.name is None:
            print("Printing info about vibrational object")
        else:
            print("Printing info on {}".format(self.name))
        print("Vibrational type: {}".format(self.type))
        print("Vibrational dimension: {}".format(self.dim))
        print("Reduced mass: {}".format(self.mu))


class HarmonicVibration(Vibration):
    """
    Class for handling molecular vibrations in the harmonic approximation
    """

    def __init__(self, Molpara, dim=3, Vname=None) -> None:
        """
        Generator routine for handling harmonic vibrations in a molecule.
        """
        try:
            super().__init__(Molpara, dim, Vtype="Harmonic", Vname=Vname)
            self.check_molpara(Molpara, ["omega", "dDdx"])
            self.omega = Molpara['omega']
            self.dDdx  = Molpara['dDdx']
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError: {}".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))

    def set_Harmonic_matrix(self) -> None:
        """
        """
        from qutip import qeye, create, destroy
        pref = self.omega
        self.H0 = pref * (create(self.dim) * destroy(self.dim) + 0.5 * qeye(self.dim))

    def set_Harmonic_transmatrix(self) -> None:
        """
        """
        from qutip import create, destroy
        pref = self.dDdx * math.sqrt(2.*self.mu*self.omega)**(-1.)
        self.Hmat = pref * (create(self.dim) + destroy(self.dim))  

    def print_info(self):
        """
        Prints info on harmonic vibrational state to screen
        """
        Vibration.print_info(self)
        print("harmonic angular frequency (omega): {}".format(self.omega))
        print("dDdX: {}".format(self.dDdx))



