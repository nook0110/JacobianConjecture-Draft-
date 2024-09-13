from enum import Enum
from aenum import Enum, NoAlias
from typing import Any
import sympy
import random
import logging
import itertools
from sympy import plot_implicit, symbols, Eq
from phcpy.dimension import get_core_count
from phcpy.solver import solve
from phcpy.solutions import filter_real
from datetime import datetime
import mpmath

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class Function:
    def __init__(self, first : sympy.poly, second : sympy.poly):
        self.first_ = first
        self.second_ = second

    def __str__(self) -> str:
        return f"x -> {self.first_}, y -> {self.second_}"
    
    def solve(self, point = [0,0]):
        H = [f'{self.first_ - point[0]};',
        f'{self.second_ - point[1]};' ]
        nbcores = get_core_count()
        solutions = []
        try:
            for solution in solve(H, tasks=1, dictionary_output=True):
                solutions.append((solution['x'], solution['y']))
        except:
            return None
        return list(set(solutions))
    
    def get_jacobian(self) -> sympy.Matrix:
        F = sympy.Matrix([self.first_, self.second_])
        return F.jacobian([x, y])
    
    def get_degree(self) -> int:
        return self.degree_

    first_ : sympy.poly
    second_ : sympy.poly
    contracts_ : bool
class CheckResult(Enum):
        _settings_ = NoAlias
        NOT_HOLDS = False
        HOLDS = True
        NON_GENERAL_POSITION = True
        CONTRACTS = True

class Checker:
    def __init__(self, function : Function) -> None:
        self.function_ = function
        self.jacobian_ = function.get_jacobian()
        self.det_ = self.jacobian_.det()
    
    def calculate(self, solutions):
        values = []
        for solution in solutions:
            val = self.det_
            val = val.subs(x, solution[0])
            val = val.subs(y, solution[1])
            values.append(val.evalf(1000))
        from itertools import combinations
        ans = sum(map(lambda x : 1/x, values))
        ans = sum((map(lambda x : sympy.N(sympy.prod(x), 10000000), combinations(values, len(values) - 1))))
        return ans 

    def calculate_for_point(self, point = [0,0]):
        solutions = self.function_.solve(point)
        print(*sorted(list(map(lambda x : (x[0], x[1]), solutions)), key = lambda x : (x[0].real, x[1].real)), sep="\n")
        return self.calculate(solutions) 
    
    function_ : Function
    jacobian_ : sympy.Matrix
    det_ : Any

x, y=sympy.symbols('x y')
from sympy.parsing.sympy_parser import parse_expr
phi_x = -10*x**2 - 3*x*y**2 + y**5
phi_y = -x**5 - y**5 - 9*y**4
degree = 25
point = (-98, -35)
print(Checker(function=Function(phi_x, phi_y)).calculate_for_point(point).evalf())
