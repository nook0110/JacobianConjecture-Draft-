from enum import Enum
from aenum import Enum, NoAlias
from typing import Any
import sympy
import random
import logging
import itertools
import main
from sympy import solve
from sympy import Eq
from sympy import evalf
x,y=sympy.symbols('x y')
phi_x = -x*y + x
phi_y = y**2 + y
point = (61, -11)
print(sympy.nsolve([Eq(phi_x, point[0]), Eq(phi_y, point[1])], [x, y], (-21 + sympy.I, 20 + sympy.I)))

def generate_random_values(size) -> list:
    return [random.randrange(-100, 100) for _ in range(size)]

print(generate_random_values(100))