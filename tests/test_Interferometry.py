#!/usr/bin/env python3
import os
import sys
sys.path.append("/home/martin/Work/Qutip/modules/")
import unittest 
import Interferometry as In
import math
    
class TestMethods(unittest.TestCase): 
    # negative test function to test if values are almost equal with place 
    omega_max = 10.
    domega    = 0.05.
    tau_max   = 20.
    dtau      = 2.
    name = "TestMethods"
    def test_Efield_interferometry_impact_2level(self): 
        decimalPlace = 12
        pass

    def test_get_dtau(self):
        dtau = In.get_dtau(self.omega_max)
        decimalPlace = 12

        self.assertAlmostEqual(dtau, 0.5*self.omega_max**(-1.), decimalPlace) 

    def test_get_taumax(self):
        decimalPlace = 12
        taumax = In.get_taumax(self.domega)
        self.assertAlmostEqual(taumax, 0.5*self.domega**(-1.), decimalPlace) 

    def test_get_ntau(self):
        ntau = In.get_ntau(self.tau_max, self.dtau)
        self.assertEqual(ntau, 11)
        
    def test_taus_from_mol(self):
        pass



