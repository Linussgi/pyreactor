import numpy as np

class Component:
    def __init__(self, name, order, sho_coeffs: list[float]):
        self.order = order
        self.name = name
        self.cp_coeffs = sho_coeffs


    def enthalpy_change(self, initial_temp, final_temp) -> int:
        t_1 = initial_temp / 1000
        t_2 = final_temp / 1000

        A, B, C, D, E = self.cp_coeffs

        term_1 = (A * t_1) + (B * t_1 ** 2 / 2) + (C * t_1 ** 3 / 3) + (D * t_1 ** 4 / 4) - (E / t_1)
        term_2 = (A * t_2) + (B * t_2 ** 2 / 2) + (C * t_2 ** 3 / 3) + (D * t_2 ** 4 / 4) - (E / t_2)

        return term_2 - term_1
    

    def entropy_change(self, initial_temp, final_temp) -> int:
        t_1 = initial_temp / 1000
        t_2 = final_temp / 1000

        A, B, C, D, E = self.cp_coeffs

        term_1 = (A * np.log(t_1)) + (B * t_1) + (C * t_1 ** 2 / 2) + (D * t_1 ** 3 / 3) - (E / (2 *  t_1 ** 2))
        term_2 = (A * np.log(t_2)) + (B * t_2) + (C * t_2 ** 2 / 2) + (D * t_2 ** 3 / 3) - (E / (2 *  t_2 ** 2))

        return term_2 - term_1