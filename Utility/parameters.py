"""
Module containing parameters used by the library
"""
import math

eps0    = (4.*math.pi)**(-1.) # vacuum permittivitiy (atomic units)
eps0_SI = 8.8541878128  #vacuum permittivitiy (SI units) ref: wikipedia
c       = 137.035999084 # speed of light (atomic units), ref: wikipedia
c_SI    = 299792458. # speed of light (SI units), ref: wikipedia
fs2au   = 41.341373 # converting femtoseconds to atomic time
mu_p2e  = 1836.152673426 # proton to electron mass ratio

Hartree2J = 4.3597447e-18 # Hartree to Joule
Bohr2m    = 5.2917721e-11 # Bohr to meter
atime2s   = 2.4188843e-17 # atomic ime to second

I0_SI     = Hartree2J / ((100.*Bohr2m)**2. * atime2s) # I0 to SI units  (I0 has SI units W / cm^2 = J / (s * cm^2))

Molecule_keys  = ['B', 'D', 'Da', 'Mol']
Molecule_types = ['MgH+', 'CaH+', 'CaD+', 'Custom']


def check_molpara(para:dict, *args):
    """
    Chek consitency of the para file
    """
    if not args:
        raise TypeError("At least one argument must be passed to compare to para")
    else:
        for arg in args:
            for i in range(len(arg)):
                if arg[i] not in Molecule_keys:
                    raise ValueError("Key must be in the list of acceptable keys {}, got {}".format(Molecule_keys, arg[i]))
                if arg[i] not in para:
                    raise KeyError("key {} not in para.".format(arg[i]))
                elif not isinstance(para[arg[i]], float):
                    if arg[i] != "Mol":
                        raise TypeError("If not Mol, {} needs to be of type float, got {}.".format(arg[i], type(para[arg[i]])))
                elif para[arg[i]] < 0.:
                    raise ValueError("{} needs to be a positive number, got {}".format(key, para[arg[i]]))


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


def kappa_from_B_sigma(B, sigma):
    """
    """
    kappa = 3. * B * 2. * math.sqrt(2.*math.log(2.)) * sigma

    return kappa
    


"""
Systems of molecular ions
"""

class Molecule():
    """
    A class containing molecular parameters
    """
    parameter_declaration = {'B': 'Rotational constant', 'Da': 'Polarizability anisotropy', 'D': 'Dipole moment', 'mu': 'Reduced mass', 'omega': 'Harmonic vibrational frequency', 're': 'Equilibrium distance', 'dDdx': 'Derivative of dipole moment w.r.t. the internuclear distance', 'a': 'First centrifugal constant', 'Mol': 'Molecular species'}


    def __init__(self, Molecule:dict|str, custom_para=None, vibration=False) -> None:
        """
        The constructor medthod
        """
        try:
            if isinstance(Molecule, str):
                if (Molecule == "Custom") and (custom_para is None):
                    raise Exception("If Custom molecule is used to initialize we must provide a dictionary with parameters, none given.")
                else:
                    Molpara = set_MolParams(Molecule, custom_para=custom_para)
            else:
                Molpara = Molecule
            if ('B' in Molpara) and ('Da' in Molpara):
                for key, val in Molpara.items():
                    if (key == 'a')  and not vibration:
                        setattr(self, key, 0.)
                    else:
                        setattr(self, key, val)


            else:
                raise KeyError("A molecule must at least have a definite rotational constant and polarizability anisotropy. Molecule has B: {}, Molecule has Da: {}".format('B' in Molpara, 'Da' in Molpara))
        except KeyError as ke:
            print("KeyError: {}".format(ke))
        except TypeError as te:
            print("TypeError: {}".format(te))
        except ValueError as ve:
            print("ValueError: {}".format(ve))
        except Exception as e:
            print("Exception: {}".format(e))

    def print_info(self, full_info=False):
       """
       Prints out the molecular parameters to screen.
       """
       for attr, val in self.__dict__.items():
           if full_info:
               if attr in Molecule.parameter_declaration:
                   print("{}: {} ({})".format(attr, val, Molecule.parameter_declaration[attr]))
               else:
                   print("{}: {}".format(attr, val))
           else:
               print("{}: {}".format(attr, val))


def set_MolParams(Molecule, custom_para=None, vibration=False) -> dict:
    """
    Retruns a dictionary of molecular paramters
    """
    if vibration:
        print("WARNING, THE VIBRATIONAL PARAMETERS HAVE NOT BEEN PROPERLY IMPLEMENTED YET")
    if Molecule == 'MgH+':
        Mol = {'B': 2.88*10.**(-5.), 'Da': 16.20, 'D': 1.18, 'a': 1.*10.**(-7.), 'Mol': 'MgH+'}
        if vibration:
            Mol.update({'mu': 4., 'omega': 0.37, 'dDdx': 1.5, 'a': 1.*10.**(-7.), 're': 2.,'Mol': 'MgH+'})
    elif Molecule == 'CaH+':
        Mol  = {'B': 2.15*10.**(-5.), 'Da': 16.59, 'D': 2.38, 'Mol': 'CaH+'}
        if vibration:
            Mol.update({'mu': 4., 'omega': 0.37, 'dDdx': 1.5, 'a': 1.*10.**(-7.), 're': 2.,'Mol': 'MgH+'})
    elif Molecule == 'CaD+':
        Mol  = {'B': 1.10*10.**(-5.), 'Da': 16.59, 'D': 2.38, 'Mol': 'CaD+'}
        if vibration:
            Mol.update({'mu': 4., 'omega': 0.37, 'dDdx': 1.5, 'a': 1.*10.**(-7.), 're': 2.,'Mol': 'MgH+'})
    elif Molecule == 'Custom':
        if custom_para is None:
            raise Exception("If Custom molecule is used to initialize we must provide a dictionary with parameters, none given.")
        else:
            Mol = {}
            for key, val in custom_para.items():
                Mol.update({key: val})
            Mol.update({'Mol': 'Custom'})
    else:
        print("Warning molecular system {} not found in list, returning empty dictionary".format(Molecule))
        Mol = {}

    return Mol

def add_vibParams(Molecule, vib_model="Harmonic"):
    """
    Add vibrational parameters to the molecule
    """
    print("Warning, vibrational parameters are not correctly implemented yet!!")
    if vib_model == "Harmonic":
        if Molecule['Mol'] == 'MgH+':
            Molecule.update({'mu': 24.*25./ (24. + 25.) * mu_p2e, 'omega': 0.5, 'dDdx': 1.15})
        elif Molecule['Mol'] == 'CaH+':
            Molecule.update({'mu': 40.*41./ (40. + 41.) * mu_p2e, 'omega': 0.5, 'dDdx': 1.15})
        elif Molecule['Mol'] == 'CaD+':
            Molecule.update({'mu': 40.*42./ (40. + 42.) * mu_p2e, 'omega': 0.5, 'dDdx': 1.15})
        else:
            raise TypeError("add_vibParams does not accept molecule {} at the moment".format(Molecule['Mol']))



def update_MolParams(Molecule, **kwargs):
    """
    Updates the dictionary Molecule with the values in kwargs
    """
    Mol = Molecule.update(kwargs)

    return Mol

def set_state(state, dim, mrep=False):
    """
    Args:
    --------
        state: int, Qobj
            The wanted molecular state

        dim: int
            Dimension of the basis (number of j-states)

        mrep: Boolean
            If true include non-zero m-states in the representation

        Returns:
        -------
            istate: Qobj
                The desired state in the desired representation
    """
    from qutip import basis

    try:
        if isinstance(state, int):
            if mrep:
                if state < (dim + 1)**2:
                    istate = basis((dim + 1)**2,state)
                else:
                    raise ValueError("Can't set to a state larger than the basis size")
            else:
                if state < dim:
                    istate = basis(dim, state)
                else:
                    raise Exception("Can't set to a state larger than the basis size")

        elif hasattr(state.type, 'ket'):
            istate = state
        else:
            print("Warning, no initial state was provided, setting to ground state!")
            if mrep:
                istate = basis((dim + 1)**2,0)
            else:
                istate = basis(dim, 0)
        return istate

    except ValueError as ve:
        print(ve)
    except Exception as e:
        print(e)



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


