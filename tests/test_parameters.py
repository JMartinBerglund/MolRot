#!/usr/bin/env python3
import os
import sys
sys.path.append("/home/martin/Work/Qutip/modules/")
import unittest 
import parameters as pm
import math
    
class TestParamters(unittest.TestCase): 
    # negative test function to test if values are almost equal with place 
    def test_eps0(self): 
        decimalPlace = 12
        # assert function() to check if values are almost equal 
        self.assertAlmostEqual(pm.eps0, (4.*math.pi)**(-1.), decimalPlace) 

    def test_c(self): 
        decimalPlace = 10
        # assert function() to check if values are almost equal 
        self.assertAlmostEqual(pm.c, 137.035999084, decimalPlace) 

    def test_fs2au(self): 
        decimalPlace = 10
        # assert function() to check if values are almost equal 
        self.assertAlmostEqual(pm.fs2au, 41.341373, decimalPlace) 


if __name__ == '__main__': 
        unittest.main() 
