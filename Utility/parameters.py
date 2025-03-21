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



"""
Pulse parameters
"""

def P_from_I0_sigma(I0, Dalpha, sigma):
    """
    Returns the integrated pulse strength from the maximum pulse intenstiy, polarizability anisotropy and    sigma
    """
    P = math.sqrt(0.5*math.pi) * Dalpha * eps0**(-1.) * c**(-1.) * I0 * sigma
    return P


def I0_from_P_sigma(P, Dalpha, sigma):
    """
    Returns the maximum pulse intensity from the integrated pulse strength, the polarizability anisotropy    and sigma
    """
    I0 = P * eps0 * c / (math.sqrt(0.5*math.pi) * Dalpha * sigma)
    return I0


def sigma_from_P_I0(P, Dalpha, I0):
    """
    Returns sigma from the maximum pulse intensity, the integrated pulse strength and the polarizability     anisotropy
    """
    sigma = P * eps0 * c / (math.sqrt(0.5*math.pi) * Dalpha * I0)
    return sigma


"""
Systems of molecular ions
"""
def set_MolParams(Molecule, **kwargs):
    """
    Retruns a dictionary of molecular paramters
    """
    if Molecule == 'MgH+':
        Mol = {'B': 2.88*10**(-5.), 'Da': 16.20, 'D': 1.18, 'Mol': 'MgH+'}
    elif Molecule == 'CaH+':
        Mol  = {'B': 2.15*10**(-5.), 'Da': 16.59, 'D': 2.38, 'Mol': 'CaH+'}
    elif Molecule == 'CaD+':
        Mol  = {'B': 1.10*10**(-5.), 'Da': 16.59, 'D': 2.38, 'Mol': 'CaD+'}
    elif Molecule == 'Custom':
        Mol = kwargs
        Mol.update({'Mol': 'Custom'})
    else:
        print("Warning molecular system {} not found in list, returning empty dictionary".format(Molecule))
        Mol = {}

    return Mol

def update_MolParams(Molecule, **kwargs):
    """
    Updates the dictionary Molecule with the values in kwargs
    """
    Mol = Molecule.update(kwargs)

    return Mol

"""
Combinations of molecular parameters and electric field
"""

def get_chiD(Molecule: float, eps0: float) -> float:
    """
    Returns the parameters chi_D from rotational constant, dipole moment and the electric field
    """
    chiD = Molecule['D'] * eps0 * Molecule['B']**(-1.)
    return chiD

def get_eps0(Molecule: dict, chiD: float) -> float:
    """
    Returns the electric field strength from rotational constant, dipole moment and the parameters chiD
    """
    eps0 = Molecule['B'] * chiD * Molecule['D']**(-1.)
    return eps0


"""
Printing information of the molecular parameters
"""

def print_molecular_info(Molecule):
    """
    Prints out what molecular system we are using
    """
    try:
        print("Molecular ion used for the optimization: {}".format(Molecule['Mol']))
    except Molecule['Mol'] == None:
        print("No molecular ion system provided")


def print_Brot(Molecule):
    """
    Prints out the rotational constant of the molecule
    """
    try:
        print("Rotational constant used for the optimization: {} (atomic units)".format(Molecule['B']))
    except Molecule == None:
        print("No molecular ion system provided")
    except Molecule['B'] not in Molecule:
        print("No rotational constant provided")

def print_Dalpha(Molecule):
    """
    Prints out the polarizability anisotropy of the molecule
    """
    try:
        print("Polarizability anisotropy used for the optimization: {} (atomic units)".format(Molecule['Da']))
    except Molecule == None:
        print("No molecular ion system provided")
    except Molecule['Da'] not in Molecule:
        print("No polarizability anisotropy provided")

def print_dip(Molecule):
    """
    Prints out the dipole moment of the molecule
    """
    try:
        print("Dipole moment used for the optimization: {} (atomic units)".format(Molecule['D']))
    except Molecule == None:
        print("No molecular ion system provided")
    except Molecule['D'] not in Molecule:
        print("No dipole moment provided")


