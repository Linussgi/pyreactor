import numpy as np
from scipy.optimize import root_scalar

from component import Component
from constants import R


class Reaction:
    def __init__(self, reactants: list[Component], products: list[Component], rxn_std_dh: float, rxn_std_ds: float, eos: str, reac_temp: float):
        self.reactants = reactants
        self.products = products
        self.rxn_std_dh = rxn_std_dh
        self.eos = eos
        self.rxn_std_ds = rxn_std_ds
        self.reac_temp = reac_temp


    def calculate_rxn_gibbs(self, *, std_temp=298):
        rxn_enthalpy = 0
        rxn_entropy = 0

        for comp in self.reactants:
            rxn_enthalpy += comp.enthalpy_change(self.reac_temp, std_temp) * comp.order
            rxn_entropy += comp.entropy_change(self.reac_temp, std_temp) * comp.order

        for comp in self.products:
            rxn_enthalpy += comp.enthalpy_change(std_temp, self.reac_temp) * comp.order
            rxn_entropy += comp.entropy_change(std_temp, self.reac_temp) * comp.order

        rxn_enthalpy += self.rxn_std_dh
        rxn_entropy += self.rxn_std_ds

        return rxn_enthalpy - self.reac_temp * rxn_entropy


    def calculate_rxn_k(self, rxn_gibbs):
        return np.exp(-rxn_gibbs / (R * self.reac_temp))
    

    def reaction_equation(self, chi, pressure, reac_init, prod_init):
        numerator = 1
        denominator = 1
        
        equil_prod = 0
        equil_reac = 0

        delta_mols = 0

        for index, comp in enumerate(self.products):
            prod_term = prod_init[index] + comp.order * chi

            numerator *= prod_term ** (comp.order)
            equil_prod += prod_term

            delta_mols += comp.order
        
        for index, comp in enumerate(self.reactants):
            reac_term = reac_init[index] - comp.order * chi

            denominator *= reac_term ** (comp.order)
            equil_reac += reac_term

            delta_mols -= comp.order

        total_mol_term = ((equil_prod + equil_reac) ** (-delta_mols))
        pressure_term = (pressure ** delta_mols)

        return (numerator / denominator) * total_mol_term  * pressure_term


    def calculate_conversion(self, pressure, reac_init, prod_init, K):
        epsilon = 1e-5
        chi_min = -max(n0 / comp.order for comp, n0 in zip(self.products, prod_init)) - epsilon
        chi_max = min(n0 / comp.order for comp, n0 in zip(self.reactants, reac_init)) + epsilon

        def residual(chi):
            try:
                val = self.reaction_equation(chi, pressure, reac_init, prod_init) - K
                if not np.isfinite(val):
                    return 1e10
                return val
            except Exception:
                return 1e10

        # Check signs at the edges:
        f_min = residual(chi_min)
        f_max = residual(chi_max)

        if f_min * f_max > 0:
            raise RuntimeError(f"No sign change in bracket: f({chi_min})={f_min}, f({chi_max})={f_max}")

        sol = root_scalar(residual, method="brentq", bracket=[chi_min, chi_max])
        if sol.converged:
            return sol.root
        else:
            raise RuntimeError("Root finding did not converge.")



