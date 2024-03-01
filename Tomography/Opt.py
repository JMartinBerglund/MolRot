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
        #print(S0, S3)
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
        #print(S0, S3)
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
    def min_DP_omega_S3(x, Ufull, Op, dim, l0, l3):
        from Tomomod import PauliN
        from Utility import UBWO, Proj
        P1 = x * (Ufull.Pulses.P[0] + Ufull.Pulses.P[1])
        P2 = (1. - x) * (Ufull.Pulses.P[0] + Ufull.Pulses.P[1])
        Ufull.update_full_operators(P=[P1, P2], Pind=[0,1])
        Op_BW = UBWO(Ufull.U, Op)
        S0 = Proj(Op_BW, PauliN(0, dim))
        S3 = Proj(Op_BW, PauliN(3, dim))

        return l0 * (S0 - 1.)**2. + l3 * S3**2.



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
                Pulses = Ut.Pulses([par['P1'], par['P2']], [par['t1'], var['t2']])
                U = Ut.EvolutionOperators(Pulses, par['B'], dim)
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
        from scipy.optimize import minimize
        
        if self.target['Target'] == 'S3':
            if self.method['Method'] == 'SP':
                opt = minimize(self.min_SP, x0=self.varlist['P1'], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
            elif self.method['Method'] == 'DP_t2':
                opt = minimize(self.min_DP_t2_S3, x0=[self.varlist['t2']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
            elif self.method['Method'] == 'DP_omega':
                opt = minimize(self.min_DP_omega_S3, x0=[self.varlist['omega']], args=(U, self.O, self.dim, self.lagrange['l0'], self.lagrange['l3']))
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
                

    def run_optimization(self) -> None:
        """
        Checks that method and target are in supported lists and calls the optimization routing
        """
        if (self.method['Method'] in self.Method_list) and (self.target['Target'] in self.Target_list):
            self.optimization()
        else:
            raise Exception("Optimization was requested with either unsuported method or target")



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


