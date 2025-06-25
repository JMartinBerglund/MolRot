"""
Module containing classes for pulses. We consider pulses in the impact approximation and without the impact approximation. also pulses for tomography and interferometry.
"""

import os
import sys
sys.path.append("/home/martin/work/qutip/modules")


class Optimizer():
    """Optimizer routine for backwards propagation of measurment operators
    
    Methods:
    
        min_SP :
            Minimizes the pulse strength with one pulse

        min_SP_t1:

    """

    Method_list = ['SP', 'SP_t1', 'DP_t1', 'DP_t2', 'DP_t1_t2', 'DP_omega', 'DP_omega_t2', 'DP_omega_t1_t2', 'DP']
    """
    List of supported methods for optimization
        
        Args:
            
            SP: string
                Single pulse optimization (target S3)
            
            SP_t1: string
                Single pulse optimizataion with time delay (targets: S1, S2)
    """
    
    Target_list = ['S1', 'S2', 'S3']
    """
    List of targets for optimization
        
        Args:

            S1: string
                Sets target to Pauli1

            S2: string
                Sets target to Pauli2

            S3: string
                Sets target to Pauli3
    """

    def __init__(self, O0, dim:int, variables:dict, params:dict, method:dict, target:dict, lagrange:float) -> None:
        """Init routine"""
        self.O         = O0
        self.dim       = dim
        self.varlist   = variables
        self.paramlist = params
        self.method    = method
        self.target    = target
        self.lagrange  = lagrange
        self.opt       = None
        #print(self.target['Target'], self.target['Target'] in self.Target_list)


    def optimization(self) -> None:
        """
        Checks for inconsisencies before performing an optimization 
        """
        met = self.method['Method']
        tar = self.target['Target']
        var = self.varlist
        par = self.paramlist

        if tar == 'S3':
            if met == 'SP':
                if 'P1' in var:
                    self.optimize()
                else:
                    raise Exception
            elif (met == 'DP_t2'):
                if ('P1' in par) and ('P2' in par) and ('t2' in var):
                    self.optimize()
                else:
                    raise Exception
            elif met == 'DP_omega':
                if ('P' in par) and ('t2' in par) and ('omega' in var):
                    self.optimize()
            else:
                print("At the moment we are only focusing on the SP, DP_t2 and DP_omega methods for target S3", met, "will hopefully be implemented soon")
                self.opt = None

        elif (tar == 'S1') or (tar == 'S2'):
            if met == 'SP':
                if ('t1' in var) and ('P1' in par):
                    self.optimize()
                else:
                    raise Exception
            elif met == 'DP_t1':
                if ('t1' in var) and ('P1' in par) and ('P2' in par) and ('t2' in par) and ('B' in par):
                    self.optimize()
            else:
                print("At the moment we are only focusing on the SP and DP_t1 methods.", met, "will hopefully be implemented soon")
                self.opt = None


    def run_optimization(self) -> None:
        """
        Checks that method and target are in supported lists and calls the optimization routing
        """
        if (self.method['Method'] in self.Method_list) and (self.target['Target'] in self.Target_list):
            self.optimization()
        else:
            raise Exception("Optimization was requested with either unsuported method or target")




class ImpactOptimizer(Optimizer):

    def __init__(self, O0, dim:int, variables:dict, params:dict, method:dict, target:dict, lagrange:float):
        super().__init__(O0, dim, variables, params, method, target, lagrange)

    @staticmethod
    def min_SP(x:float, Up, Op, dim:int, l0:float, l3:float) -> float:
        """Minimizes the projection of $\sigma_3$ for backwards propagation with a single pulse, with the restraint that the total population of the 2-level system remains 1
        
        Args:

            x : float
                The pulse fluence to optimize
            
            Up : QuTiP Qobj
                Evolution operator of the pulse in the impulse approximation
            
            Op : QuTiP Qobj
                The measurment operator to be propagated
            
            dim : int
                dimension of the basis $\ge 2$
            
            l0 : float
                Lagrange multiplyer for projection on $\sigma_0$
            
            l3 : float
                Lagrange multiplyer for projection on $\sigma_3$

        Outputs:
            
            L : float
                Lagrange function to minimize
        """
        import Utility as Ut
        from Utility import UBWO, Proj
        from Tomomod import PauliN

        Up.update_pulse_operator(x)
        Op_BW = UBWO(Up.Up, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S3 = Proj(Op_BW, PauliN(3, dim))
        L = l0 * (S0 - 1.)**2. + l3 * S3**2.

        return L


    @staticmethod
    def min_SP_S1(x:float, U, Op, dim:int, l0:float, l1:float) -> float:
        """Minimizes the projection of :math:´\\sigma_1´ for backwards propagation with a single pulse
        
        Args:
            
            x : float
                The pulse fluence to optimize
            
            Up : QuTiP Qobj
                Evolution operator of the pulse in the impulse approximation
            
            Op : QuTiP Qobj
                The measurment operator to be propagated
            
            dim : int
                dimension of the basis $\ge 2$
            
            l0 : float
                Lagrange multiplyer for projection on $\sigma_0$
            
            l1 : float
                Lagrange multiplyer for projection on $\sigma_1$

        Returns:
            
            L : float
                Lagrange function to minimize
        """
        import Utility as Ut
        from Utility import UBWO, Proj
        from Tomomod import PauliN
        U.update_full_operator(which="t", value=x[0])
        Op_BW = UBWO(U.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S1 = Proj(Op_BW, PauliN(1, dim))
        L = l0 * (S0 - 1.)**2. + l1 * (S1 - 1.)**2.
 

        return L

    @staticmethod
    def min_SP_S2(x:float, U, Op, dim:int, l0:float, l2:float) -> float:
        """Minimizes the projection of $\sigma_2$ for backwards propagation with a single pulse
        
        Args:
            
            x : float
                The pulse fluence to optimize
            
            Up : QuTiP Qobj
                Evolution operator of the pulse in the impulse approximation
            
            Op : QuTiP Qobj
                The measurment operator to be propagated
            
            dim : int
                dimension of the basis $\ge 2$
            
            l0 : float
                Lagrange multiplyer for projection on $\sigma_0$
            
            l1 : float
                Lagrange multiplyer for projection on $\sigma_1$

        Returns:
            
            L : float
                Lagrange function to minimize

        """
        import Utility as Ut
        from Utility import UBWO, Proj
        from Tomomod import PauliN
        U.update_full_operator(which="t", value=x[0])
        Op_BW = UBWO(U.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S2 = Proj(Op_BW, PauliN(2, dim))
        #print(S0, S3)
        L = l0 * (S0 - 1.)**2. + l2 * (S2 - 1.)**2.

        return L



    @staticmethod
    def min_DP_t2_S1(x:float, Ufull, Op, dim:int, l0:float, l1:float) -> float:
        """Minimizes the projection of $\sigma_1$ for backwards propagation with two pulses with a time delay.
        
        Args:
            
            x : float
                The pulse fluence to optimize
            
            Ufull : QuTiP Qobj
                Evolution operator of the pulse in the impulse approximation
            
            Op : QuTiP Qobj
                The measurment operator to be propagated
            
            dim : int
                dimension of the basis $\ge 2$
            
            l0 : float
                Lagrange multiplyer for projection on $\sigma_0$
            
            l1 : float
                Lagrange multiplyer for projection on $\sigma_1$

        Returns:
            
            L : float
                Lagrange function to minimize

        """
        import Utility as Ut
        from Tomomod import UBWO, Proj, PauliN
        Ufull.update_full_operators(time=x, tind=[1])
        Op_BW = UBWO(Ufull.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S1 = Proj(Op_BW, PauliN(1, dim))
        #print('In min_DP_t2_S1')
        #print(x)
        L = l0 * (S0 - 1.)**2. + l1 * (S1 - 1.)**2.

        return L

    @staticmethod
    def min_DP_t2_S2(x, Ufull, Op, dim, l0, l2):
        import Utility as Ut
        from Tomomod import UBWO, Proj, PauliN
        Ufull.update_full_operators(time=x, tind=[1])
        Op_BW = UBWO(Ufull.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S2 = Proj(Op_BW, PauliN(2, dim))
        #print('In min_DP_t2_S2')
        #print(x)

        return l0 * (S0 - 1.)**2. + l2 * (S2 - 1.)**2.


    @staticmethod
    def min_DP_t1_S1(x, Ufull, Op, dim, l0, l1):
        import Utility as Ut
        from Utility import UBWO, Proj
        from Tomomod import PauliN
        Ufull.update_full_operators(time=x, tind=[0])
        Op_BW = UBWO(Ufull.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S1 = Proj(Op_BW, PauliN(1, dim))

        return l0 * (S0 - 1.)**2. + l1 * (S1 - 1.)**2.

    @staticmethod
    def min_DP_t1_S2(x, Ufull, Op, dim, l0, l2):
        import Utility as Ut
        from Tomomod import PauliN
        from Utility import UBWO, Proj
        Ufull.update_full_operators(time=x, tind=[0])
        Op_BW = UBWO(Ufull.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S2 = Proj(Op_BW, PauliN(2, dim))

        return l0 * (S0 - 1.)**2. + l2 * (S2 - 1.)**2.


    @staticmethod
    def min_DP_t2_S3(x, Ufull, Op, dim, l0, l3):
        from Tomomod import PauliN
        from Utility import UBWO, Proj
        Ufull.update_full_operators(time=x, tind=[1])
        Op_BW = UBWO(Ufull.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S3 = Proj(Op_BW, PauliN(3, dim))

        return l0 * (S0 - 1.)**2. + l3 * S3**2.

    @staticmethod
    def min_DP_omega_S3(x, Ufull, P, Op, dim, l0, l3):
        from Tomomod import PauliN
        from Utility import UBWO, Proj
        P1 = x * P
        P2 = (1. - x) * P
        Ufull.update_full_operators(P=[P1, P2], Pind=[0,1])
        Op_BW = UBWO(Ufull.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S1 = Proj(Op_BW, PauliN(1, dim))
        S2 = Proj(Op_BW, PauliN(2, dim))
        S3 = Proj(Op_BW, PauliN(3, dim))
        
        # make sure that omega is in [0, 1]
        #if (omega < 0.) or (omega > 1.):
        #    return 100.
        #else:
        #    return l0 * (S0 - 1.)**2. + l3 * S3**2.
        #lag = l0 * (S0 - 1.)**2. + l3 * S3**2.  
        #lag = l0 * S0**2. + (1. -(S1**2. + S2**2.))**2. + l3 * S3**2.  
        lag = l0 * S0**2. + l3 * S3**2.  
        return lag 



    @staticmethod
    def set_evolution_operator(dim:int, met:dict, tar:dict, var=None, par=None):
        """Set the evolution operator based on pulse strength(s), rotational constant and time delay(s).

            Args:
                dim: int
                    Dimension of the basis
                method: dict
                    Which optimization method to use
                tar: dict
                    Target of the optimization
                var: dictionary
                    The variables of the optimization routine
                par: dict
                    The parameters of the optimization routine

            Returns:
                U : Qobj
                    The evolution operator
        """
        import Utility as Ut
        if tar == 'S3':
            if met == 'SP':
                U = Ut.ImpulseEvolutionOperator(var['P1'], dim) 
                return U
            elif met == 'DP_t2':
                Pulsepara = {'Ps': [par['P1'], par['P2']], 'taus': [par['t1'], var['t2']]}
                #Pulses = Ut.Pulses([par['P1'], par['P2']], [par['t1'], var['t2']])
                Pulses = Ut.ImpactPulses(Pulsepara)
                U = Ut.EvolutionOperators(Pulses, par['B'], dim)
                return U
            elif met == 'DP_omega':
                P1 = var['omega'] * par['P']
                P2 = (1. - var['omega']) * par['P']
                Pulsepara = {'Ps': [P1, P2], 'taus': [par['t1'], par['t2']]}
                Pulses = Ut.ImpactPulses(Pulsepara)
                U = Ut.EvolutionOperators(Pulses, par['B'], dim)
                return U
            else:
                raise RuntimeError("WARNING:", met, 'not implemented yet for S3!')

        elif (tar == 'S1') or (tar == 'S2'):
            if met == 'SP':
                U = Ut.FullEvolutionOperator(par['P1'], par['B'], var['t1'], dim) 
                return U
            elif met == 'DP_t1':
                Pulses = Ut.Pulses([par['P1'], par['P2']], [var['t1'], par['t2']])
                U = Ut.EvolutionOperators(Pulses, par['B'], dim)
                return U
            else:
                raise RuntimeError("WARNING:", met, 'not implemented yet for S1, S2!')
        else:
            raise RuntimeError("WARNING:", met, 'not implemented yet for S1, S2!')


    def minimizer(self, U):
        """Method to select optimization method based on the target and method dictionaries

            Args:
                U: Qobj
                    Evolution operator

            Returns:
                opt: Intance of result from minimize
        """
        from scipy.optimize import minimize, Bounds
        
        if self.target['Target'] == 'S3':
            if self.method['Method'] == 'SP':
                opt = minimize(self.min_SP, x0=self.varlist['P1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
            elif self.method['Method'] == 'DP_t2':
                opt = minimize(self.min_DP_t2_S3, x0=[self.varlist['t2']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
            elif self.method['Method'] == 'DP_omega':
                opt = minimize(self.min_DP_omega_S3, x0=[self.varlist['omega']], args=(U, self.paramlist['P'], self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']), bounds=Bounds(lb=0., ub=1.))
            else:
                raise RuntimeError("No valid method given for S3 minimization:", self.method['Method'])
        
        elif self.target['Target'] == 'S1':
            if self.method['Method'] == 'SP':
                opt = minimize(self.min_SP_S1, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
            elif self.method['Method'] == 'DP_t1':
                opt = minimize(self.min_DP_t1_S1, x0=[self.varlist['t1']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
            else:
                raise RuntimeError("No valid method given for S1 minimization:", self.method['Method'])
        
        elif self.target['Target'] == 'S2':
            if self.method['Method'] == 'SP':
                opt = minimize(self.min_SP_S2, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l2']))
            elif self.method['Method'] == 'DP_t1':
                opt = minimize(self.min_DP_t1_S2, x0=[self.varlist['t1']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l2']))
            else:
                raise RuntimeError("No valid method given for S2 minimization:", self.method['Method'])

        else:
            raise Exception

        return opt
                

    def optimize(self) -> None:
        """Sets the evolution operator base on target and method and calls the minimizer to execute the optimization
        """
        from scipy.optimize import minimize
        import Utility as Ut

        U = self.set_evolution_operator(self.dim, self.method['Method'], self.target['Target'], self.varlist, self.paramlist)
        opt = self.minimizer(U)
        self.opt = opt
        #if self.method['Method'] == 'SP':
        #    if self.target['Target'] == 'S3':
        #        U = Ut.ImpulseEvolutionOperator(self.varlist['P1'], self.dim)
        #        opt = minimize(self.min_SP, x0=self.varlist['P1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
        #        self.opt = opt
        #    elif self.target['Target'] == 'S1':
        #        U = Ut.FullEvolutionOperator(self.paramlist['P1'], self.paramlist['B'], self.varlist['t1'], self.dim)
        #        opt = minimize(self.min_SP_S1, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
        #        self.opt = opt
        #    elif self.target['Target'] == 'S2':
        #        U = Ut.FullEvolutionOperator(self.paramlist['P1'], self.paramlist['B'], self.varlist['t1'], self.dim)
        #        opt = minimize(self.min_SP_S2, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l2']))
        #        self.opt = opt
        #elif self.method['Method'] == 'SP_t1':
        #    U = Ut.FullEvolutionOperator(self.paramlist['P1'], self.B, self.varlist['t1'], self.dim)
        #    if self.target == 'S1':
        #        opt = minimize(self.min_SP_t1, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
        #    elif self.target == 'S2':
        #        opt = minimize(self.min_SP_t1, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
#
        #    self.opt = opt

        #elif self.method['Method'] == 'DP_t2':
        #    if self.target['Target'] == 'S1':
        #        Pulses = Ut.Pulses([self.paramlist['P1'], self.paramlist['P2']], [0., self.varlist['t2']])
        #        U = Ut.EvolutionOperators(Pulses, self.paramlist['B'], self.dim)
        #        opt = minimize(self.min_DP_t2_S1, x0=[self.varlist['t2']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
        #        self.opt = opt
        #    elif self.target['Target'] == 'S2':
        #        Pulses = Ut.Pulses([self.paramlist['P1'], self.paramlist['P2']], [0., self.varlist['t2']])
        #        U = Ut.EvolutionOperators(Pulses, self.paramlist['B'], self.dim)
        #        opt = minimize(self.min_DP_t2_S2, x0=[self.varlist['t2']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l2']))
        #        self.opt = opt
        #    elif self.target['Target'] == 'S3':
        #        print(self.varlist['t2'])
        #        print(self.paramlist['P1'], self.paramlist['P2'], self.lagrange['l0'], self.lagrange['l3'])
                #print(self._min_DP_t2_S3(self.varlist['t2'], self.paramlist['P1'], self.paramlist['P2'], self.lagrange['l0'], self.lagrange['l3']))
        #        self.test(self.varlist['t2'])
        #        opt, tol = minimize(self._min_DP_t2_S3, x0=[self.varlist['t2']], args=(self.paramlist['P1'], self.paramlist['P2'], self.lagrange['l0'], self.lagrange['l3']))
                #opt, tol = minimize(self.test, x0=3.)
        #    else:
        #        raise Exception('Unrecognized target:', self.target['Target'], 'for optimization', self.method['Method'])
        #else:
        #    raise Exception('Method:', self.method['Method'], 'not implemented yet')


        #rint('In optimize')

    #def optimize_SP_S3(self):
    #    print('In optimize_SP_S3')

    #def optimize_DP_t2_S3(self):
    #    print('In optimize_DPt2_S3')

        #if (met == 'SP'):
        #    if 'P1' in var:
        #        if tar == 'S3':
        #            self.optimize()
        #        else:
        #            raise Exception(tar,'not compatible with variable pulse strength for method', met)
        #    elif 't1' in var:
        #        if (tar == 'S1') or (tar == 'S2'):
        #            self.optimize()
        #        else:
        #            raise Exception(tar,'not compatible with variable time for method', met)

        #    else:
        #        raise Exception('One of the required variables P1 or t1 not provided to single pulse optimization')
        #elif met == 'DP_t2':
        #    if 't2' in var:
        #        if ('P1' in par) and ('P2' in par):
        #            if (tar == 'S1') or (tar == 'S2') or (tar == 'S3'):
        #                self.optimize()
        #            else:
        #                raise Exception(tar,'not yet implemented for', met)
        #        else:
        #            raise Exception('Required parameters P1 and P2 not provided to double pulse optimization with variable t2')
        #    else:
        #        raise Exception('Required variable t2 not provided to double pulse optimization with variable t2')
        #else:
        #    raise Exception('Implementation work in progress. No supported optimization method was provided')
                


class GaussOptimizer(Optimizer):

    def __init__(self, O0, dim:int, variables:dict, params:dict, method:dict, target:dict, lagrange:float, options=None):
        """
        """
        super().__init__(O0, dim, variables, params, method, target, lagrange)
        self.options=options # Options for the sesolver. Might be necessarry for convergence


    @staticmethod
    def min_SP(x:float, H0, HI, Da, sigma, tmax, tmin, dt, Op, dim:int, l0:float, l3:float, options) -> float:
        """Minimizes the projection of $\sigma_3$ for backwards propagation with a single pulse, with the restraint that the total population of the 2-level system remains 1
        
        Args:

            x : float
                The pulse fluence to optimize
            
            Up : QuTiP Qobj
                Evolution operator of the pulse in the impulse approximation
            
            Op : QuTiP Qobj
                The measurment operator to be propagated
            
            dim : int
                dimension of the basis $\ge 2$
            
            l0 : float
                Lagrange multiplyer for projection on $\sigma_0$
            
            l3 : float
                Lagrange multiplyer for projection on $\sigma_3$

        Outputs:
            
            L : float
                Lagrange function to minimize
        """
        #mport Utility as Ut
        from parameters import I0_from_P_sigma as I0ps
        from Utility import Gauss_me, Proj
        from Tomomod import PauliN
        from qutip import sesolve, basis

        # Update the pulse intensity (Equivalent to the pulse fluence at constant sigma) Up.update_pulse_operator(x)
        I0 = I0ps(x, Da, sigma)
        istate1 = basis(dim, 0)
        istate2 = basis(dim, 1)
        res1 = sesolve([H0, [HI, Gauss_me]], istate1, np.arange(tmax, tmin, -dt), args={}, options=options)
        #Op_BW = UBWO(Up.Up, Op)
        #S0 = Proj(Op_BW, PauliN(0, dim))
        #S3 = Proj(Op_BW, PauliN(3, dim))
        L = 1. #l0 * (S0 - 1.)**2. + l3 * S3**2.
        print(H0)
        print(HI)

        return L

 


    def optimize(self):
        """
        """
        from scipy.optimize import minimize
        import Utility as Ut

        #U = self.set_evolution_operator(self.dim, self.method['Method'], self.target['Target'], self.varlist, self.paramlist)
        self.opt = self.minimizer()
        #self.opt = opt



    def minimizer(self):
        """Method to select optimization method based on the target and method dictionaries

            Args:
                U: Qobj
                    Evolution operator

            Returns:
                opt: Intance of result from minimize
        """
        from scipy.optimize import minimize, Bounds
        
        if self.target['Target'] == 'S3':
            if self.method['Method'] == 'SP':
                opt = minimize(self.min_SP, x0=self.varlist['P1'], args=(self.H0, self.Hint, self.Molecule.Da, self.Pulses.sigma, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
            elif self.method['Method'] == 'DP_t2':
                opt = minimize(self.min_DP_t2_S3, x0=[self.varlist['t2']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
            elif self.method['Method'] == 'DP_omega':
                opt = minimize(self.min_DP_omega_S3, x0=[self.varlist['omega']], args=(U, self.paramlist['P'], self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']), bounds=Bounds(lb=0., ub=1.))
            else:
                raise RuntimeError("No valid method given for S3 minimization:", self.method['Method'])
        
        elif self.target['Target'] == 'S1':
            if self.method['Method'] == 'SP':
                opt = minimize(self.min_SP_S1, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
            elif self.method['Method'] == 'DP_t1':
                opt = minimize(self.min_DP_t1_S1, x0=[self.varlist['t1']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l1']))
            else:
                raise RuntimeError("No valid method given for S1 minimization:", self.method['Method'])
        
        elif self.target['Target'] == 'S2':
            if self.method['Method'] == 'SP':
                opt = minimize(self.min_SP_S2, x0=self.varlist['t1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l2']))
            elif self.method['Method'] == 'DP_t1':
                opt = minimize(self.min_DP_t1_S2, x0=[self.varlist['t1']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l2']))
            else:
                raise RuntimeError("No valid method given for S2 minimization:", self.method['Method'])

        else:
            raise Exception

        return opt
 



class Optimizer_BCH(Optimizer):
    """
    Optimization class in the Baker-Campbell-Hausdorff correction to the impulse approximation. 
    """

    def __init__(self, O0, dim:int, variables:dict, params:dict, method:dict, target:dict, lagrange:float) -> None:
        """
        """
        super().__init__(O0, dim, variables, params, method, target, lagrange)

    @staticmethod
    def set_evolution_operator(dim:int, met:dict, tar:dict, var=None, par=None):
        """Set the evolution operator based on pulse strength(s), rotational constant and time delay(s).

            Args:
                dim: int
                    Dimension of the basis
                method: dict
                    Which optimization method to use
                tar: dict
                    Target of the optimization
                var: dictionary
                    The variables of the optimization routine
                par: dict
                    The parameters of the optimization routine

            Returns:
                U : Qobj
                    The evolution operator
        """
        import Utility as Ut
        #from parameteres import set_kappa
        from Utility import U2U1, set_BCH_corr
        if tar == 'S3':
            if met == 'SP':
                Pulsepara = {'Ps': [var['P1']], 'taus': [0.], 'sigma': par['sigma']}
                Pulses = Ut.Pulses(Pulsepara, Ptype="impulse_BCH")
                U = Ut.EvolutionOperators_CBH(Pulses, par['B'], 2)
                #U  = U2U1(set_BCH_corr(par['B'], var['P1'], par['sigma']), UI.Up)
                #print("Tried correction")
                return U
            elif met == 'DP_t2':
                Pulses = Ut.Pulses([par['P1'], par['P2']], [par['t1'], var['t2']])
                UI = Ut.EvolutionOperators(Pulses, par['B'], dim)
                U  = U2U1(set_BCH_corr(par['B'], var['P1'], par['sigma']), UI)
                return U
            elif met == 'DP_omega':
                P1 = var['omega'] * par['P']
                P2 = (1. - var['omega']) * par['P']
                Pulses = Ut.Pulses([P1, P2], [par['t1'], par['t2']])
                U = Ut.EvolutionOperators(Pulses, par['B'], dim)
                return U
            else:
                raise RuntimeError("WARNING:", met, 'not implemented yet for S3!')

        elif (tar == 'S1') or (tar == 'S2'):
            if met == 'SP':
                UI = Ut.FullEvolutionOperator(par['P1'], par['B'], var['t1'], dim) 
                U  = U2U1(set_BCH_corr(par['B'], var['P1'], par['sigma']), UI)
                return U
            elif met == 'DP_t1':
                Pulses = Ut.Pulses([par['P1'], par['P2']], [var['t1'], par['t2']])
                UI = Ut.EvolutionOperators(Pulses, par['B'], dim)
                U  = U2U1(set_BCH_corr(par['B'], var['P1'], par['sigma']), UI)
                return U
            else:
                raise RuntimeError("WARNING:", met, 'not implemented yet for S1, S2!')
        else:
            raise RuntimeError("WARNING:", met, 'not implemented yet for S1, S2!')

    @staticmethod
    def min_SP(x:float, Up, Op, dim:int, l0:float, l3:float) -> float:
        """Minimizes the projection of $\sigma_3$ for backwards propagation with a single pulse, with the restraint that the total population of the 2-level system remains 1
        
        Args:

            x : float
                The pulse fluence to optimize
            
            Up : QuTiP Qobj
                Evolution operator of the pulse in the impulse approximation
            
            Op : QuTiP Qobj
                The measurment operator to be propagated
            
            dim : int
                dimension of the basis $\ge 2$
            
            l0 : float
                Lagrange multiplyer for projection on $\sigma_0$
            
            l3 : float
                Lagrange multiplyer for projection on $\sigma_3$

        Outputs:
            
            L : float
                Lagrange function to minimize
        """
        import Utility as Ut
        from Utility import UBWO, Proj
        from Tomomod import PauliN

        Up.update_pulse_operator(x)
        Op_BW = UBWO(Up.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S3 = Proj(Op_BW, PauliN(3, dim))
        L = l0 * (S0 - 1.)**2. + l3 * S3**2.

        return L


class Optimizer2P(Optimizer):
    """Optimizer wrapper for a 2 pulse evolution operator"""

    # INPUTS
    # @U0:         Initial guess for the evolution operator (EvolutioOperator instance)
    # @O0:         Measurment operators to be used in the optimization
    # @variables:  Independent variables to be varied by the optimizer(dict)
    # @methods:    Which Lagrangian to use, and value of the Lagrangian multipliers (dict)
    def __init__(self, U0, O0, variables, params, methods):
        """Initializing an optimizer"""
        self.U          = U0
        self.Ops        = O0
        self.varlist    = variables
        self.paramlist  = params
        self.method     = methods['Method']
        self.l0         = methods['Lag0']
        self.l1         = methods['Lag1']
        self.l2         = methods['Lag2']
        self.l3         = methods['Lag3']


    def optimize(self):
        from scipy.optimize import minimize

        #result = minimize(self.method, )

        pass


class OptimizerEfield():
    """
    Optimizer for measuring the electric field within a Paul trap
    """

    def __init__(self, Mol, method, approximation):
        """ Initialziation of the optimizer

        Args:
        -----

        Mol: string
            The molecular ion we are using for the optimization

        method: string
            Which methdod we are using for optimization
        """
        import parameters as pm

    
        self.Mol    = pm.set_MolParams(Mol)
        self.method = method
        self.approx = approximation


    def print_mol_info(self):
        """
        """
        import parameters as pm
        pm.print_molecular_info(self.Mol)


    def print_B(self):
        """
        """
        import parameters as pm
        pm.print_Brot(self.Mol)

    def print_Da(self):
        """
        """
        import parameters as pm
        pm.print_Dalpha(self.Mol)

    def print_D(self):
        """
        """
        import parameters as pm
        pm.print_dip(self.Mol)

    def print_full_mol_info(self):
        """
        """
        print("")
        self.print_mol_info()
        self.print_B()
        self.print_Da()
        self.print_D()
        print("")


    def print_method_info(self):
        """
        """
        try:
            print("Method used for the optimization: {}".format(self.method))
        except self.method == None:
            print("No method provided")

    def print_approximation_info(self):
        """
        """
        try:
            print("Approximation used for the optimization: {}".format(self.approx))
        except self.approx == None:
            print("No approximation provided")

    def print_tau_info(self):
        """
        """
        try:
            print("Tau-grid used for the optimization: taumin: {}, taumax: {}, dtau: {}, ntau: {}".format(self.tauv[0], self.tauv[-1], self.tauv[1] - self.tauv[0], len(self.tauv)))
        except self.tauv == None:
            print("No tau-grid provided")



    def print_P_info(self):
        """
        """
        try:
            print("Pulse fluence used for the optimization: {}".format(self.P))
        except self.P == None:
            print("No pulse fluence provided")

    def print_eps0_ref_info(self):
        """
        """
        try:
            print("Efield to be found by the optimization: {}".format(self.eps0_ref))
        except self.eps0_ref == None:
            print("No reference field provided")

class OptimizerEfield_ls(OptimizerEfield):
    """
    """

    def __init__(self, Mol, approximation, tauv, P, reference_spectrum, reference_omega, eps0_guess=0., eps0_ref=None, n=None):
        """
        """
        super().__init__(Mol, method="least_squares", approximation=approximation)
        self.refspec = reference_spectrum
        self.refom   = reference_omega
        self.tauv    = tauv
        self.P       = P
        self.eps0    = eps0_guess
        self.ep0_ref = eps0_ref
        self.n       = n

    

    @staticmethod
    def opt_imp_2l(x, B, D, P, tau, reference_spectrum, reference_omega):
        """
        """
        from Interferometry import Efield_interferometry_impact_2level as EII2L 
        from numpy.fft import fft, fftfreq
        import numpy as np
        from Utility import least_square_fitting as ls
        eps0 = x[0]
        ltau = len(tau)
        dtau = tau[1] - tau[0]

        # Obtain the interferogram in the 2level impact approximation
        c_sq = EII2L(B, eps0, D, P, tau)
        
        # Get the spectrum
        Fc_sq = fft(c_sq) 
        spec = np.abs(Fc_sq[0:ltau//2]) / ltau
        omega = fftfreq(ltau, dtau)[0:ltau//2]


        return ls(spec, reference_spectrum)
    
    @staticmethod
    def opt_imp(x, B, D, P, tau, reference_spectrum, reference_omega, n):
        """
        """
        from Interferometry import ImpactEfieldInterferometry as IEI 
        from numpy.fft import fft, fftfreq
        import numpy as np
        from Utility import least_square_fitting as ls
        from Utility import get_DipEigen as DE
        eps0 = x[0]
        print("Testing opt_imp, printing current eps0:", eps0)
        
        # Making sure that the electric field is positive
        if eps0 < 0.:
            return 10000.
        else:
            # Get the initial state from diagonalization of the dipole Hamiltonian
            _, Eigstates = DE(B, D, eps0, n) # HERE WE NEED TO CONTUNIE!!!
            istate = Eigstates[0]

            # Obtain the interferogram in the impact approximation
            IntImpEf = IEI(Ps=[P, P], ts=tau, B=B, eps=eps0, D=D, istate=istate, dim=n, Name="")
            IntImpEf.run_interferometry()

            # Get the spectrum from interferogram
            IntImpEf.get_spectrum()
            return ls(IntImpEf.spectrum, reference_spectrum) 
 


    def optimize_Efield(self, fullInfo=False):
        """
        """
        from scipy.optimize import minimize as mini
        from scipy.optimize import Bounds as bd

        print("")
        print("Running a least square fit to reference spectrum to obtain the electric field")
        print("Approximation method used is {}".format(self.approx))

        if self.approx == "Impact_2l":
            res = mini(self.opt_imp_2l, x0=self.eps0, args=(self.Mol['B'], self.Mol['D'], self.P, self.tauv, self.refspec, self.refom))
            if fullInfo:
                return res
            else:
                return res.x
        elif self.approx == "Impact_full":
            if self.n is not None:
                if True: # change to optinal 2-level approximation
                    # Runs a 2-level optimization first to get a reasonable initial value for eps0 for the numerical optimization
                    print("Testing optimize_Efield, method Impact_full. Entering 2-level appr.")
                    res = mini(self.opt_imp_2l, x0=self.eps0, args=(self.Mol['B'], self.Mol['D'], self.P, self.tauv, self.refspec, self.refom))
                    eps0_2l = res.x[0]
                    print("Testing optimize_Efield, method Impact_full. Entering full basis appr.")
                    if eps0_2l == 0.:
                        uppb = 0.0001
                    else:
                        uppb = 10.*eps0_2l
                    lowb = -uppb

                    res = mini(self.opt_imp, x0=eps0_2l, args=(self.Mol['B'], self.Mol['D'], self.P, self.tauv, self.refspec, self.refom, self.n), bounds=bd(lb=lowb, ub=uppb), options={"maxiter": 8, 'disp': True})
                    print("Testing optimize_Efield, method Impact_full. Exiting full basis appr.")
                    if fullInfo:
                        return res
                    else:
                        return res.x
                else:
                    res = mini(self.opt_imp, x0=self.eps0, args=(self.Mol['B'], self.Mol['D'], self.P, self.tauv, seld.dim, self.refspec, self.refom))
                    if fullInfo:
                        return res
                    else:
                        return res.x
            else:
                print("Need a basis size for the full impact approximation!, None given.")
        else:
            print("Approximation: {} not yet implemented!".format(self.approx))
            

