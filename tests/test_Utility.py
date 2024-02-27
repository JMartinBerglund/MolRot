#!/usr/bin/env python3
import os
import sys
sys.path.append("/home/martin/Work/Qutip/modules/")
import unittest 
import Utility as Ut
import math
    
class TestPulses(unittest.TestCase): 
    # negative test function to test if values are almost equal with place 
    P = [1., 2.]
    t = [3., 4.]
    name = "TestPulses"
    def test_init(self): 
        decimalPlace = 12
        Pulses = Ut.Pulses(P=self.P, t=self.t, name=self.name)
        # assert function() to check if values are almost equal 
        for i in range(len(self.P)):
            self.assertAlmostEqual(Pulses.P[i], self.P[i], decimalPlace) 
            self.assertAlmostEqual(Pulses.t[i], self.t[i], decimalPlace) 

        self.assertEqual(Pulses.Ptype, "impulse") 
        self.assertEqual(Pulses.nP, len(self.P)) 
        self.assertEqual(Pulses.name, self.name) 

    def test_copy_pulse(self): 
        Pulses = Ut.Pulses(P=self.P, t=self.t, name=self.name)
        Pulses2 = Pulses.copy_pulse()
        decimalPlace = 10
        for i in range(len(self.P)):
            self.assertAlmostEqual(Pulses.P[i], Pulses2.P[i], decimalPlace) 
            self.assertAlmostEqual(Pulses.t[i], Pulses2.t[i], decimalPlace) 

        self.assertEqual(Pulses.Ptype, Pulses2.Ptype) 
        self.assertEqual(Pulses.nP, Pulses2.nP) 
        self.assertEqual(Pulses.name, Pulses2.name) 

class TestEvolutionOperator(unittest.TestCase): 
    # negative test function to test if values are almost equal with place 
    dim    = 2
    Otype  = "impulse"
    Otype2 = "impulse2"
    name   = "TestOperator"
    name2  = "TestOperator2"
    odd    = False
    def test_EO(self): 
        decimalPlace = 12
        # Test init
        Eo = Ut.EvolutionOperator(self.dim, self.Otype, self.name, self.odd)
        self.assertEqual(Eo.dim, self.dim)
        self.assertEqual(Eo.Otype, self.Otype)
        self.assertEqual(Eo.name, self.name)
        self.assertEqual(Eo.odd, self.odd)

        # Test rename
        Eo.rename_operator(self.name2)
        self.assertEqual(Eo.name, self.name2)

        # Test update
        Eo.update_operator_type(self.Otype2)
        self.assertEqual(Eo.Otype, self.Otype2)

class TestImpulseEvolutionOperator(unittest.TestCase): 
    # negative test function to test if values are almost equal with place 
    dim    = 2
    P      = 3.
    P2     = 1.
    Otype  = "pulse"
    name   = "TestImpulseOperator"
    odd    = False
    def test_IEO(self): 
        decimalPlace = 12
        # Test init
        IEo = Ut.ImpulseEvolutionOperator(self.P, self.dim, self.name, self.odd)
        UP  = Ut.UI(self.P, self.dim, self.odd)
        UP2 = Ut.UI(self.P2, self.dim, self.odd)
        Project = Ut.Proj(UP-IEo.Up, UP-IEo.Up)
        self.assertEqual(IEo.P, self.P)
        self.assertEqual(IEo.dim, self.dim)
        self.assertEqual(IEo.Otype, self.Otype)
        self.assertEqual(IEo.name, self.name)
        self.assertEqual(IEo.odd, self.odd)
        self.assertAlmostEqual(Project, 0., decimalPlace)

        # Test updates
        IEo.update_pulse(self.P2)
        self.assertAlmostEqual(IEo.P, self.P2, decimalPlace)
        IEo.update_pulse_operator(self.P2)
        Project2 = Ut.Proj(UP2-IEo.Up, UP2-IEo.Up)
        self.assertAlmostEqual(Project2, 0., decimalPlace)

class TestFreeEvolutionOperator(unittest.TestCase): 
    # negative test function to test if values are almost equal with place 
    dim    = 2
    B      = 3.
    B2     = 1.
    t      = 0.5
    t2     = 1.
    Otype  = "free"
    name   = "TestFreeOperator"
    odd    = False
    which  = 'B'
    which2 = 't'
    def test_FEO(self): 
        decimalPlace = 12
        # Test init
        FEo = Ut.FreeEvolutionOperator(self.B, self.t, self.dim, self.name, self.odd)
        UF  = Ut.Uf(self.B, self.t, self.dim, self.odd)
        UF2 = Ut.Uf(self.B, self.t2, self.dim, self.odd)
        Project = Ut.Proj(UF-FEo.Uf, UF-FEo.Uf)
        self.assertEqual(FEo.B, self.B)
        self.assertEqual(FEo.t, self.t)
        self.assertEqual(FEo.dim, self.dim)
        self.assertEqual(FEo.Otype, self.Otype)
        self.assertEqual(FEo.name, self.name)
        self.assertEqual(FEo.odd, self.odd)
        self.assertAlmostEqual(Project, 0., decimalPlace)

        # Test updates
        FEo.update_B(self.B2)
        self.assertAlmostEqual(FEo.B, self.B2, decimalPlace)
        FEo.update_time(self.t2)
        self.assertAlmostEqual(FEo.t, self.t2, decimalPlace)
        FEo.update_free_operator(self.which, self.B)
        Project2 = Ut.Proj(UF2-FEo.Uf, UF2-FEo.Uf)
        self.assertAlmostEqual(Project2, 0., decimalPlace)
        FEo.update_free_operator(self.which2, self.t)
        Project3 = Ut.Proj(UF-FEo.Uf, UF-FEo.Uf)
        self.assertAlmostEqual(Project3, 0., decimalPlace)

class TestHamiltonians(unittest.TestCase):

    from qutip import basis
    dim   = 2
    dim2  = 3
    B     = 1.
    odd   = False
    Htype = 'oper'
    s1    = basis(dim, 0)
    s2    = basis(dim, 1)
    t1    = basis(dim2, 0)
    t2    = basis(dim2, 1)
    t3    = basis(dim2, 2)
    Da    = 2.

    def test_H0s(self):
        decimalPlace = 12
        
        H1 = Ut.H0(self.B, self.dim, self.odd)
        self.assertEqual(H1.shape[0], 2)
        self.assertEqual(H1.shape[1], 2)
        self.assertEqual(H1.type, self.Htype)
        self.assertAlmostEqual(H1.matrix_element(self.s1, self.s1), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H1.matrix_element(self.s1, self.s2), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H1.matrix_element(self.s2, self.s1), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H1.matrix_element(self.s2, self.s2), 6.*self.B+0.j, decimalPlace)

        H2 = Ut.H0(self.B, self.dim, True)
        self.assertEqual(H2.shape[0], 2)
        self.assertEqual(H2.shape[1], 2)
        self.assertEqual(H2.type, self.Htype)
        self.assertAlmostEqual(H2.matrix_element(self.s1, self.s1), 2.*self.B+0.j, decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s1, self.s2), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s2, self.s1), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s2, self.s2), 12.*self.B+0.j, decimalPlace)

        H3 = Ut.H0(self.B, self.dim2, self.odd)
        self.assertEqual(H3.shape[0], 3)
        self.assertEqual(H3.shape[1], 3)
        self.assertEqual(H3.type, self.Htype)
        self.assertAlmostEqual(H3.matrix_element(self.t1, self.t1), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t1, self.t2), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t1, self.t3), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t2, self.t1), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t2, self.t2), 6.*self.B+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t2, self.t3), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t3, self.t1), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t3, self.t2), 0.+0.j, decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.t3, self.t3), 20.*self.B+0.j, decimalPlace)


    def test_HIs(self):
        from parameters import eps0, c
        decimalPlace = 12
        
        H1 = Ut.get_HIntMat(self.dim2, self.odd)
        self.assertEqual(H1.shape[0], 3)
        self.assertEqual(H1.shape[1], 3)
        self.assertAlmostEqual(H1[0,0], 3.**(-1.), decimalPlace)
        self.assertAlmostEqual(H1[1,1], 11./21., decimalPlace)
        self.assertAlmostEqual(H1[2,2], 39./77., decimalPlace)
        self.assertAlmostEqual(H1[0,1], 2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H1[1,0], 2./(3.*math.sqrt(5.)), decimalPlace)

        H2 = Ut.H1(self.Da, self.dim, self.odd)
        self.assertEqual(H2.shape[0], 2)
        self.assertEqual(H2.shape[1], 2)
        self.assertAlmostEqual(H2.matrix_element(self.s1,self.s1), -0.25*self.Da/3., decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s1,self.s2), -0.25*self.Da*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s2,self.s1), -0.25*self.Da*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s2,self.s2), -0.25*self.Da*11./21., decimalPlace)

        H3 = Ut.H1_I0(self.Da, self.dim, self.odd)
        self.assertEqual(H3.shape[0], 2)
        self.assertEqual(H3.shape[1], 2)
        self.assertAlmostEqual(H3.matrix_element(self.s1,self.s1), -self.Da/(2.*eps0*c)/3., decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.s1,self.s2), -self.Da/(2.*eps0*c)*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.s2,self.s1), -self.Da/(2.*eps0*c)*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.s2,self.s2), -self.Da/(2.*eps0*c)*11./21., decimalPlace)

    def test_HImpact(self):
        decimalPlace = 12
        P1 = 0.
        P2 = 2.

        H1 = Ut.HI(P1, self.dim, self.odd)
        H2 = Ut.HI(P2, self.dim, self.odd)
        H3 = Ut.HImi(P2, self.dim, self.odd)
        self.assertEqual(H1.shape[0], 2)
        self.assertEqual(H1.shape[1], 2)
        self.assertAlmostEqual(H1.matrix_element(self.s1,self.s1), 0., decimalPlace)
        self.assertAlmostEqual(H1.matrix_element(self.s1,self.s2), 0., decimalPlace)
        self.assertAlmostEqual(H1.matrix_element(self.s2,self.s1), 0., decimalPlace)
        self.assertAlmostEqual(H1.matrix_element(self.s2,self.s2), 0., decimalPlace)
        
        #Needs more thought
        self.assertEqual(H2.shape[0], 2)
        self.assertEqual(H2.shape[1], 2)
        self.assertAlmostEqual(H2.matrix_element(self.s1,self.s1), -P2/3., decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s1,self.s2), -P2*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s2,self.s1), -P2*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H2.matrix_element(self.s2,self.s2), -P2*11./21., decimalPlace)
        
        self.assertEqual(H3.shape[0], 2)
        self.assertEqual(H3.shape[1], 2)
        self.assertAlmostEqual(H3.matrix_element(self.s1,self.s1), 1.j*P2/3., decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.s1,self.s2), 1.j*P2*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.s2,self.s1), 1.j*P2*2./(3.*math.sqrt(5.)), decimalPlace)
        self.assertAlmostEqual(H3.matrix_element(self.s2,self.s2), 1.j*P2*11./21., decimalPlace)
        

class testUs(unittest.TestCase):
    from qutip import basis
    dim   = 2
    dim2  = 3
    B     = 1.
    odd   = False
    Utype = 'oper'
    s1    = basis(dim, 0)
    s2    = basis(dim, 1)
    ID2   = s1*s1.dag() + s2*s2.dag()
    t1    = basis(dim2, 0)
    t2    = basis(dim2, 1)
    t3    = basis(dim2, 2)
    ID3   = t1*t1.dag() + t2*t2.dag() +t3*t3.dag()
    Da    = 2.

    def test_UI(self):
        import numpy as np
        import cmath
        from scipy.linalg import inv
        from qutip import Qobj
        decimalPlace = 12
        P1 = 0.
        P2 = 2.
        a  = 1./3.
        b  = 2./(3.*math.sqrt(5.))
        d  = 11./21.
        summ = a+d
        diff = a-d
        nuplus = 0.5*(summ + math.sqrt(diff**2. + 4.*b**2.))
        numinus = 0.5*(summ - math.sqrt(diff**2. + 4.*b**2.))
        lplus  = P2 * nuplus
        lminus = P2 * numinus
        aplus = a - nuplus
        nplus  = math.sqrt(aplus**2. + b**2.)

        V = nplus**(-1.) * np.array([[aplus, -b], [b, aplus]])
        L = np.array([[cmath.exp(1.j*lminus), 0.], [0., cmath.exp(1.j*lplus)]])
        Vinv = inv(V)
        Mat = np.matmul(V, np.matmul(L,Vinv))
        U3_comp = Qobj(Mat)


        U1 = Ut.UI(P1, self.dim, self.odd)
        U2 = Ut.UI(P1, self.dim2, self.odd)
        U3 = Ut.UI(P2, self.dim, self.odd)
        
        # Assert U1
        self.assertEqual(U1.shape[0], 2)
        self.assertEqual(U1.shape[1], 2)
        self.assertEqual(U1.type, self.Utype)
        diff = U1 - self.ID2
        self.assertAlmostEqual(Ut.Proj(diff, diff), 0., decimalPlace)

        # Assert U2
        self.assertEqual(U2.shape[0], 3)
        self.assertEqual(U2.shape[1], 3)
        self.assertEqual(U2.type, self.Utype)
        diff2 = U2 - self.ID3
        self.assertAlmostEqual(Ut.Proj(diff2, diff2), 0., decimalPlace)

        # Assert U3
        self.assertEqual(U3.shape[0], 2)
        self.assertEqual(U3.shape[1], 2)
        self.assertEqual(U3.type, self.Utype)
        diff3 = U3 - U3_comp
        self.assertAlmostEqual(Ut.Proj(diff3, diff3), 0., decimalPlace)

    def test_Uf(self):
        import numpy as np
        import cmath
        from scipy.linalg import inv
        from qutip import Qobj
        decimalPlace = 12
        B1 = 2.
        t1 = 3.

        Mat1 = np.array([[1., 0.], [0., cmath.exp(-1.j*6.*B1*t1)]])
        Mat2 = np.array([[1., 0., 0.], [0., cmath.exp(-1.j*6.*B1*t1), 0.], [0., 0., cmath.exp(-1.j*20.*B1*t1)]])
        U1 = Ut.Uf(B1, t1, self.dim, self.odd)
        U2 = Ut.Uf(B1, t1, self.dim2, self.odd)
        self.assertEqual(U1.shape[0], 2)
        self.assertEqual(U1.shape[1], 2)
        self.assertEqual(U2.shape[0], 3)
        self.assertEqual(U2.shape[1], 3)

        U1_comp = Qobj(Mat1)
        U2_comp = Qobj(Mat2)
        diff1 = U1 - U1_comp
        self.assertAlmostEqual(Ut.Proj(diff1, diff1), 0., decimalPlace)
        diff2 = U2 - U2_comp
        self.assertAlmostEqual(Ut.Proj(diff2, diff2), 0., decimalPlace)


    def test_U2U1(self):
        import pytest
        import numpy as np
        import cmath
        from scipy.linalg import inv
        from qutip import Qobj
        decimalPlace = 12
        B1 = 2.
        t1 = 3.
        P1 = 1.
        U1 = Ut.Uf(B1, t1, self.dim, self.odd)
        U2 = Ut.UI(P1, self.dim, self.odd)
        U3 = Ut.Uf(B1, t1, self.dim2, self.odd)

        U = Ut.U2U1(U1, U2)
        U_comp = U1*U2
        diff = U - U_comp
        self.assertAlmostEqual(Ut.Proj(diff, diff), 0., decimalPlace)

    def test_UFW(self):
        from qutip import Qobj, basis
        decimalPlace = 12
        P1 = 1.
        s1 = basis(self.dim, 0)
        O = s1 * s1.dag()
        U = Ut.UI(P1, self.dim, self.odd)
        A = Ut.UFW(U, s1)
        B = Ut.UFW(U, O)

        diff1 = A - B
        self.assertAlmostEqual(Ut.Proj(diff1, diff1), 0., decimalPlace)

        B_comp = U*O*U.dag()
        diff2 = B - B_comp
        self.assertAlmostEqual(Ut.Proj(diff2, diff2), 0., decimalPlace)
        

    def test_UFWO(self):
        pass

    def test_UBW(self):
        from qutip import Qobj, basis
        from Tomomod import Pauli
        decimalPlace = 12
        P1 = 1.
        Ps = 2.6726819190994378
        Pr1 = -0.31943828249997
        Pr2 = -0.947607082958686
        Pr3 = 0.

        s1 = basis(self.dim, 0)
        O = s1 * s1.dag()
        U = Ut.UI(P1, self.dim, self.odd)
        Us = Ut.UI(Ps, self.dim, self.odd)
        A = Ut.UBW(U, s1)
        B = Ut.UBW(U, O)
        C = Ut.UBW(Us, O)

        diff1 = A - B
        self.assertAlmostEqual(Ut.Proj(diff1, diff1), 0., decimalPlace)

        B_comp = U.dag()*O*U
        diff2 = B - B_comp
        self.assertAlmostEqual(Ut.Proj(diff2, diff2), 0., decimalPlace)
        
        self.assertAlmostEqual(Ut.Proj(C, Pauli(0)), 1., decimalPlace)
        self.assertAlmostEqual(Ut.Proj(C, Pauli(1)), Pr1, decimalPlace)
        self.assertAlmostEqual(Ut.Proj(C, Pauli(2)), Pr2, decimalPlace)
        self.assertAlmostEqual(Ut.Proj(C, Pauli(3)), Pr3, decimalPlace)


    def test_proj(self):
        from qutip import Qobj, basis
        from Tomomod import Pauli
        decimalPlace = 12

        # The diagonals should equal 2.
        for i in range(4):
            self.assertAlmostEqual(Ut.Proj(Pauli(i), Pauli(i)), 2., decimalPlace)
            for j in range(3):
                k = (i+j+1) % 4
                self.assertAlmostEqual(Ut.Proj(Pauli(i), Pauli(k)), 0., decimalPlace)
        




 

        #with pytest.raises(TypeError) as excinfo:  
        #    Ut.U2U1(U1, U3)  
        #assert str(excinfo.value) == "Incompatible Qobj shapes"  

 


        #self.assertAlmostEqual(H2.matrix_element(self.t1,self.t3), 0., decimalPlace)
        #self.assertEqual(H1.shape[0], 2)
        #self.assertEqual(H1.shape[1], 2)
        #self.assertEqual(H1.type, self.Htype)
        #self.assertAlmostEqual(H1.matrix_element(self.s1, self.s1), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H1.matrix_element(self.s1, self.s2), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H1.matrix_element(self.s2, self.s1), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H1.matrix_element(self.s2, self.s2), 6.*self.B+0.j, decimalPlace)

        #H2 = Ut.H0(self.B, self.dim, True)
        #self.assertEqual(H2.shape[0], 2)
        #self.assertEqual(H2.shape[1], 2)
        #self.assertEqual(H2.type, self.Htype)
        #self.assertAlmostEqual(H2.matrix_element(self.s1, self.s1), 2.*self.B+0.j, decimalPlace)
        #self.assertAlmostEqual(H2.matrix_element(self.s1, self.s2), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H2.matrix_element(self.s2, self.s1), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H2.matrix_element(self.s2, self.s2), 12.*self.B+0.j, decimalPlace)

        #H3 = Ut.H0(self.B, self.dim2, self.odd)
        #self.assertEqual(H3.shape[0], 3)
        #self.assertEqual(H3.shape[1], 3)
        #self.assertEqual(H3.type, self.Htype)
        #self.assertAlmostEqual(H3.matrix_element(self.t1, self.t1), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t1, self.t2), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t1, self.t3), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t2, self.t1), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t2, self.t2), 6.*self.B+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t2, self.t3), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t3, self.t1), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t3, self.t2), 0.+0.j, decimalPlace)
        #self.assertAlmostEqual(H3.matrix_element(self.t3, self.t3), 20.*self.B+0.j, decimalPlace)


        # Test update
        #Eo.update_operator_type(self.Otype2)
        #self.assertEqual(Eo.Otype, self.Otype2)



        # assert function() to check if values are almost equal 
        #for i in range(len(self.P)):
        #    self.assertAlmostEqual(Pulses.P[i], self.P[i], decimalPlace) 
        #    self.assertAlmostEqual(Pulses.t[i], self.t[i], decimalPlace) 

        #self.assertEqual(Pulses.Ptype, "impulse") 
        #self.assertEqual(Pulses.name, self.name) 

    #def test_copy_pulse(self): 
    #    Pulses = Ut.Pulses(P=self.P, t=self.t, name=self.name)
    #    Pulses2 = Pulses.copy_pulse()
    #    decimalPlace = 10
    #    for i in range(len(self.P)):
    #        self.assertAlmostEqual(Pulses.P[i], Pulses2.P[i], decimalPlace) 
    #        self.assertAlmostEqual(Pulses.t[i], Pulses2.t[i], decimalPlace) 

    #    self.assertEqual(Pulses.Ptype, Pulses2.Ptype) 
    #    self.assertEqual(Pulses.name, Pulses2.name) 


    #    # assert function() to check if values are almost equal 
    #    self.assertAlmostEqual(pm.c, 137.035999084, decimalPlace) 

    #def test_fs2au(self): 
    #    decimalPlace = 10
    #    # assert function() to check if values are almost equal 
    #    self.assertAlmostEqual(pm.fs2au, 41.341373, decimalPlace) 


if __name__ == '__main__': 
        unittest.main() 
