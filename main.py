from reaction import Reaction
from component import Component


ammonia = Component("NH3", 2, [0, 0, 0, 0, 0])
hydrogen = Component("H2", 3, [0, 0, 0, 0, 0])
nitrogen = Component("N2", 1, [0, 0, 0, 0, 0])

k_eq = 1e-5

haber_process = Reaction([nitrogen, hydrogen], [ammonia], 0, 0, "ideal", 600)

result = haber_process.calculate_conversion(100, [100, 300], [0], k_eq)

print(result)