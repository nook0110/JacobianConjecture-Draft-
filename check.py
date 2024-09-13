from enum import Enum
from aenum import Enum, NoAlias
from typing import Any
import sympy
import random
import logging
import itertools
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

class Function:
    def __init__(self, first : sympy.poly, second : sympy.poly, degree : int):
        self.first_ = first
        self.second_ = second
        self.degree_ = degree

    def __str__(self) -> str:
        return f"x -> {self.first_}, y -> {self.second_}"
            
    def solve(self, point = [0,0]):
        from sympy import solve
        from sympy import Eq
        
        return solve([Eq(self.first_, point[0]), Eq(self.second_, point[1])], [x, y])
    
    def get_jacobian(self) -> sympy.Matrix:
        F = sympy.Matrix([self.first_, self.second_])
        return F.jacobian([x, y])
    
    def get_degree(self) -> int:
        return self.degree_

    first_ : sympy.poly
    second_ : sympy.poly
    degree_ : int

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
    
    @staticmethod
    def is_finite_solutions(solutions) -> bool:
        from sympy import sympify
        for value in solutions:
            expr = sympify(value)
            symbols = expr.free_symbols
            if symbols:
                return False
        return True
    
    def calculate(self, solutions):
        values = []
        for solution in solutions:
            val = self.det_
            print(val)
            val = val.subs(x, solution[0].evalf())
            val = val.subs(y, solution[1].evalf())
            values.append(val.evalf())
        from itertools import combinations
        ans = sum((map(lambda x : sympy.N(sympy.prod(x)), combinations(values, len(values) - 1))))
        return ans 

    def calculate_for_point(self, point = [0,0]):
        solutions = self.function_.solve(point)
        print(*list(map(lambda x : (x[0].evalf(), x[1].evalf()), solutions)), sep="\n")
        return self.calculate(solutions) 
    
    function_ : Function
    jacobian_ : sympy.Matrix
    det_ : Any

x, y=sympy.symbols('x y')
from sympy.parsing.sympy_parser import parse_expr
phi_x = -x**2 + x
phi_y = x*y + y
degree = 2
point = (-97, -13)
print(Checker(function=Function(phi_x, phi_y, degree)).calculate_for_point(point).evalf())
