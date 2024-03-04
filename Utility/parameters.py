"""
Module containing parameters used by the library
"""
import math

eps0    = (4.*math.pi)**(-1.) # vacuum permittivitiy (atomic units)
eps0_SI = 8.8541878128  #vacuum permittivitiy (SI units) ref: wikipedia
c       = 137.035999084 # speed of light (atomic units), ref: wikipedia
c_SI    = 299792458. # speed of light (SI units), ref: wikipedia
fs2au   = 41.341373 # converting femtoseconds to atomic time

Hartree2J = 4.3597447e-18 # Hartree to Joule
Bohr2m    = 5.2917721e-11 # Bohr to meter
atime2s   = 2.4188843e-17 # atomic ime to second

I0_SI     = Hartree2J / (Bohr2m**2. * atime2s) # I0 to SI units 

