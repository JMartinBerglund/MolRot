"""Module containing molecular classes and methods with rovibrational couplings.
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
from qutip import tensor as tp




class Tensor():
    """
    """

    def __init__(self, dim1=2, dim2=2, Tname=None) -> None:
        """
        """
        try:
            self.init_dims(dim1, dim2)
            self.init_name(Tname)
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError {}:".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))

    def init_dims(self, dim1, dim2) -> None:
        """
        """
        self.dim = np.zeros(3, dtype=int)
        if isinstance(dim1, int):
            if dim1 > 0:
                self.dim[0] = dim1
            else:
                raise ValueError("The first dimension given is negative: {}, must be positive".format(dim1))
        else:
            raise TypeError("The first dimension must be a positive integer, got: {}".format(dim_1))
        if isinstance(dim2, int):
            if dim2 > 0:
                self.dim[1] = dim2
            else:
                raise ValueError("The second dimension given is negative: {}, must be positive".format(dim2))
        else:
            raise TypeError("The second dimension must be a positive integer, got: {}".format(dim2))

        self.dim[2] = dim1 * dim2

 
    def init_name(self, name):
        """
        """            
        if (name is None) or isinstance(name, str):
            self.name = name
        else:
            raise TypeError("The nema of the tensor object must be of either None or string format, got: {}",format(type(name)))



    def print_info(self) -> None:
        """
        Prints out info of the tensor object
        """

        if self.name is None:
            print("Info about nameless Tensor object")
        else:
            print("Info about Tensor object {}".format(self.name))
        print("Dimensions:")
        print("---------------------")
        print("First dimension: {}, second dimension: {}, full dimension: {}".format(self.dim[0], self.dim[1], self.dim[2]))
        print()


    def tensor_prod(self, Qobj_list):
        """
        Takes the tensor product of the Qobjs in the list

        Args:
        ------

            Qobj_list: numpy array
                List of Qobjs to apply the tensor product to

        Returns:
        -------
            tens_prod: Qobj
                The tensor product of the Qobjs in Qobj_list
            
        """
        from qutip import tensor

        tens_prod = tensor(Qobj_list)

        return tens_prod



class RovibTensor(Tensor):
    """
    """
    def __init__(self, Molpara:dict|str, dim_r=3, dim_v=2, name=None, custom_para=None, vibration=False) -> None:
        """
        """
        from parameters import Molecule as M

        try:
            super().__init__(dim_r, dim_v, name)
            self.Rot = M(Molpara, custom_para, vibration)
            self.type = "Rovib"
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError {}:".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))


    def print_info(self, full_info=False) -> None:
        """
        Prints out info of the tensor object
        """
        super().print_info()
        print("This is a {} object. It is used for rovibrational coupling.".format(self.type))
        print()
        print("Molecular info:")
        print("The molecular object has information on both rotatinoal and vibrational parameters.")
        print("Atomic untis are used throughout.")
        print("---------------------")
        self.Rot.print_info(full_info)
        print()



class ImpactHarmonicTensor(RovibTensor):
    """
    """

    def __init__(self, Molpara:dict|str, dim_r=3, dim_v=2, name=None, custom_para=None, Efield=None) -> None:
        """
        """
        try:
            if (Efield is not None) and ('D' not in Molpara):
                raise TypeError("Efield interaction requires a dipole moment!!")
            else:
                super().__init__(Molpara, dim_r, dim_v, name)
                self.approx = []
                self.approx.append("Impulse")
                self.approx.append("Harmonic")
                self.Efield = Efield
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError {}:".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))

    def set_Rotational_Hamitonians(self) -> None:
        """
        """


    def print_info(self, full_info=False) -> None:
        """
        Prints out info of the tensor object
        """
        super().print_info(full_info)
        print("Approximations:")
        print("---------------")
        print("Interaction approximation: {}".format(self.approx[0]))
        print("Vibrational approximation: {}".format(self.approx[1]))
        print()
        if self.Efield is not None:
            print("Electric field:")
            print("---------------")
            print("Field strength: {}, polar angle: {}, azimuthal angle: {}".format(self.Efield[0], self.Efield[1], self.Efield[2]))
        if hasattr(self, Rotham):
            print("Rotational Hamiltonians:")
            print("------------------------")
            print("Rotational free Hamiltonian: {}".format(self.Rotham.type()))
            print("")



class RotationTensor(Tensor):
    """
    Two coupled rotational molecules
    """

    def __init__(self, Molpara:str|dict, Molpara2=None, dim=3, istate1=0, istate2=0, name=None, custom_para=None, custom_para2=None) -> None:
        """
        Generator of the class
        """
        from parameters import Molecule as M

        try:
            super().__init__(dim, dim, name)
            self.Mol1 = M(Molpara, custom_para)
            self.H01  = None
            self.H02  = None
            self.HDa1 = None
            self.HDa2 = None
            self.Uf   = None
            self.Ui11 = None
            self.Ui12 = None
            self.Ui21 = None
            self.Ui22 = None
            if (isinstance(Molpara2, str)) or (isinstance(Molpara2, dict)):
                self.Mol2 = M(Molpara2, custom_para2)
            else:
                self.Mol2 = M(Molpara, custom_para)
            self.type = "Coupled rotation"
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError {}:".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))


    def print_info(self, full_info=False) -> None:
        """
        Prints out info of the tensor object
        """
        super().print_info()
        print("Molecular info:")
        print("Atomic untis are used throughout.")
        print("---------------------")
        print()
        print("Molecule 1:")
        self.Mol1.print_info(full_info)
        if self.Mol1.Mol != self.Mol2.Mol:
            print()
            print("Molecule 2:")
            self.Mol2.print_info(full_info)
        else:
            print("Molecule 2 is the same as molecule 1")
        print()
        if self.H01 is not None:
            print("Hamiltonian information:")
            print("------------------------")
            print()
            print("Free Hamiltonian1: {}".format(self.H01))
            print()
            print("Free Hamiltonian2: {}".format(self.H02))
            print()
        if self.Uf is not None:
            print("Impact evolution operators:")
            print("---------------")
            print()
            print("Free evolution operator: {}".format(self.Uf))
            print()
            print("Impact pulse operator1: {}".format(self.Ui11))
            print()
            print("Impact pulse operator2: {}".format(self.Ui12))
            print()
            print("Impact pulse operator3: {}".format(self.Ui21))
            print()
            print("Impact pulse operator4: {}".format(self.Ui22))
            print()
        #print("Molecule1: {}, Molecule2: {}".format(self.Mol1['Mol'], self.Mol2['Mol']))


    def set_Ham(self):
        """
        """
        from qutip import qeye
        # The identity matrix
        Id = qeye(self.dim[0])


        # First the free Hamiltonians
        H01, H02 = self.set_H0()
        #self.H01 = H01
        #self.H02 = H02
        self.H01 = self.tensor_prod([H01, Id])
        self.H02 = self.tensor_prod([Id,H02])

        # Then the polarizability interaction
        HDa1, HDa2 = self.set_HDa()
        self.HDa1 = self.tensor_prod([HDa1, Id])
        self.HDa2 = self.tensor_prod([Id, HDa2])



    def set_H0(self):
        """
        """
        from Utility import H0
        H01 = H0(self.Mol1.B, self.dim[0]) 
        H02 = H0(self.Mol2.B, self.dim[0])

        return H01, H02


    def set_HDa(self):
        """
        """
        from Utility import H1_I0 as Hint

        Hint1 = Hint(self.Mol1.Da, self.dim[0])
        Hint2 = Hint(self.Mol2.Da, self.dim[0])

        return Hint1, Hint2

    def sum_exp(self, Ham1, Ham2, const1:float, const2:float):
        O1 = -1.j * Ham1 * const1
        O2 = -1.j * Ham2 * const2
        O = O1 + O2
        U = O.expm()

        return U

    def exp(self, Ham, const:float):
        O = -1.j * Ham * const
        U = O.expm()

        return U
        


    def set_impact_U(self, tau=0., P1=0., P2=0.) -> None:
        """
        """

        if self.H01 is None:
            self.set_Ham()

        self.Uf   = self.sum_exp(self.H01, self.H02, tau, tau)
        self.Ui11 = self.exp(self.HDa1, -P1)
        self.Ui12 = self.exp(self.HDa1, -P2)
        self.Ui21 = self.exp(self.HDa2, -P1)
        self.Ui22 = self.exp(self.HDa2, -P2)

    def CX_psi(self, psi, control=2):
        """
        Args:
        -----
            psi: Qobj (ket)
                The wave function on which to operate

            control: int
                Which state is used as control, default 2
        """
        CX = self.set_CX()
        phi = CX * psi
        return phi


    def set_CX(self):
        """
        """
        from qutip import qeye, tensor, ket2dm

        Id = tensor(qeye(self.dim[0]), qeye(self.dim[0]))
        O1 = ket2dm(tensor(basis(self.dim[0],0), basis(self.dim[0],1)))
        O2 = ket2dm(tensor(basis(self.dim[0],1), basis(self.dim[0],1)))
        O3 = tensor(basis(self.dim[0],0), basis(self.dim[0],1)) * tensor(basis(self.dim[0],1), basis(self.dim[0],1)).dag()
        
        CX = Id - O1 - O2 + O3 + O3.dag()

        return CX

