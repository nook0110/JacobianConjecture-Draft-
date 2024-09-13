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

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)

x,y=sympy.symbols('x y')

def generate_random_values(size) -> list:
    return [random.randrange(-100, 100) for _ in range(size)]

def generate_random_value():
    return generate_random_values(1)[0]

def generate_random_points(size) -> list:
    return list(zip(generate_random_values(size), generate_random_values(size)))

def generate_random_point():
    return generate_random_points(1)[0]



class Function:
    def __init__(self, first : sympy.poly, second : sympy.poly):
        self.first_ = first
        self.second_ = second
        self.degree_ = 0
        self.evaluate_potential_field_extension_degree()

    def __str__(self) -> str:
        return f"x -> {self.first_}, y -> {self.second_}"
            
    def solve(self, point = [0,0]):
        
        H = [f'{self.first_ - point[0]};',
        f'{self.second_ - point[1]};' ]
        nbcores = get_core_count()
        solutions = []
        for solution in solve(H, tasks=1, dictionary_output=True):
            solutions.append((solution['x'], solution['y']))
        return list(set(solutions))
    
    def get_jacobian(self) -> sympy.Matrix:
        F = sympy.Matrix([self.first_, self.second_])
        return F.jacobian([x, y])
    
    def get_degree(self) -> int:
        return self.degree_

    def evaluate_potential_field_extension_degree(self):
        logger.info(f"Evaluating field extension for {self}")
        points_checked = 0
        max_checks = 30
        while points_checked < max_checks:
            point = generate_random_point()
            solutions = self.solve(point)
            logger.info(f"Checked {points_checked} of {max_checks}")
            points_checked += 1
            self.degree_ = max(self.degree_, len(solutions))
            if self.degree_ >= sympy.total_degree(self.first_) * sympy.total_degree(self.second_):
                print(solutions, self.degree_)
                return
    
    first_ : sympy.poly
    second_ : sympy.poly
    degree_ : int

class Generator:
    def generate(self) -> Function:
        return Function()

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
            val = val.subs(x, solution[0])
            val = val.subs(y, solution[1])
            values.append(val.evalf())
        from itertools import combinations
        ans = sum((map(lambda x : sympy.N(sympy.prod(x)), combinations(values, len(values) - 1))))
        return ans 


    def calculate_for_point(self, point = [0,0]):
        solutions = function.solve(point)
        return self.calculate(solutions) 
    
    def check_solution(self, solutions) -> bool:
        if len(solutions) > self.function_.get_degree():
            logger.error("wrong extension degree!")
        return len(solutions) >= self.function_.get_degree()
    
    def generate_possible_contraction_points(self) -> list:
        ans = []
        solutions_counter = 0
        amount = 20
        while solutions_counter < amount:
            value = generate_random_value()
            det_without_x = self.det_.subs(x, value)
            solution = sympy.solve(det_without_x, [y])
            if not solution:
                continue
            ans.append([value, solution[0]])
            
        solutions_counter = 0
        while solutions_counter < amount:
            value = generate_random_value()
            det_without_y = self.det_.subs(y, value)
            solution = sympy.solve(det_without_y, [x])
            if not solution:
                continue
            ans.append([solution[0], value])
        
        return ans
   
    def check_for_contraction(self, point) -> bool:
        solutions = self.function_.solve(point)
        return not self.is_finite_solutions(solutions)
   
 
    def test_point(self, point = [0,0]) -> CheckResult:
        logger.info(f"Checking point: {point}")
        solutions = self.function_.solve(point)
        if not self.is_finite_solutions(solutions):
            return CheckResult.CONTRACTS
        if not self.check_solution(solutions):
            logger.info(f"Amount of solutions: {len(solutions)} (expected {self.function_.get_degree()})")
            return CheckResult.NON_GENERAL_POSITION
        ev = sympy.sympify(self.calculate(solutions))
        if abs(ev.evalf()) < 0.001:
            return CheckResult.HOLDS
        logger.info(f"Test failed! Eval = {ev.evalf()}")
        return CheckResult.NOT_HOLDS 

    def test(self, iterations = 10):
        logger.info(f"Checking function: {self.function_}")
        general_position_checks = 0
        while general_position_checks < iterations:
            point = generate_random_point()
            result = self.test_point(point)
            logger.info(f"Test result for {point} is: {result.name}")
            if result == CheckResult.CONTRACTS:
                logger.info(f"Contracts at point: {point}")
                return (CheckResult.CONTRACTS, point)
            if result.value == False:
                logger.info(f"Lemma doesn't hold for function: {self.function_}")
                return (CheckResult.NOT_HOLDS, point)
            if result != CheckResult.NON_GENERAL_POSITION:
                general_position_checks += 1
                
        logger.info(f"Lemma holds for function: {self.function_}")
        return (CheckResult.HOLDS, None)
                
    function_ : Function
    jacobian_ : sympy.Matrix
    det_ : Any

        
if __name__ == '__main__':
    monomials = []
    max_deg = 2
    for deg in range(2, max_deg + 1):
        for x_deg in range(0, deg + 1):
            monomials += [x**x_deg * y**(deg - x_deg)]


    possible_monomials = list(set((itertools.combinations(monomials + [0, 0], 3))))
    coeffs_list = list(itertools.product([-1, 1], repeat = 3))

    def create(coeffs, monomials):
        return coeffs[0] * monomials[0] + coeffs[1] * monomials[1] + coeffs[2] * monomials[2]

    class CheckInfo:
        def __init__(self, function : Function, result : CheckResult, point : tuple, check_amount : int) -> None:
            self.function_ = function
            self.result_ = result
            self.point_ = point
            self.check_amount_ = check_amount

        function_ : Function
        result_ : CheckResult
        point_ : tuple
        check_amount_ : int

    import mysql.connector as con
    connection = con.connect(
          host='localhost',
          user='root',
          password='root',
          database='checked_functions'
    )
    mycursor = connection.cursor()

    def write_result(info : CheckInfo):
        sql = f'''INSERT INTO functions (Phi_x, Phi_y, ExtensionDegree, CheckResult, CheckAmount, Point) VALUES (
        \'{info.function_.first_}\', 
        \'{info.function_.second_}\',
        {info.function_.get_degree()},
        \'{info.result_.name}\',
        {info.check_amount_},
        '{info.point_}'
        );
        '''
        mycursor.execute(sql)
        connection.commit()
    
    def write_not_holds(info : CheckInfo, value : str):
        sql = f'''INSERT INTO functions (Phi_x, Phi_y, ExtensionDegree, CheckResult, CheckAmount, Point, Value) VALUES (
        \'{info.function_.first_}\', 
        \'{info.function_.second_}\',
        {info.function_.get_degree()},
        \'{info.result_.name}\',
        {info.check_amount_},
        '{info.point_}',
        '{value}'
        );
        '''
        mycursor.execute(sql)
        connection.commit()

    def check(phi_x, phi_y) -> bool:
        sql = f"SELECT EXISTS(SELECT * FROM functions WHERE Phi_x='{phi_x}' AND Phi_y='{phi_y}')"
        mycursor.execute(sql)
        return mycursor.fetchall()[0][0]

    i = 0
    to_check = []
    for coeffs_x in coeffs_list:
        for monomials_phi_x in possible_monomials:
            for coeffs_y in coeffs_list:
                for monomials_phi_y in possible_monomials:
                    i += 1
                    if i % 1000 == 0:
                        print(i/(len(coeffs_list) ** 2 * len(possible_monomials) ** 2) * 100, '%')
                    phi_x = x + create(coeffs_x, monomials_phi_x)
                    phi_y = y + create(coeffs_y, monomials_phi_y)
                    if check(phi_x, phi_y):
                        logger.info(f"Found function x->{phi_x}, y->{phi_y} in db!")
                        continue
                    to_check.append((phi_x, phi_y))
    to_check = list(set(to_check))
    random.shuffle(to_check)
    amount_checked = 0
    for phi_x, phi_y in to_check:
        logger.info(f'Left to check: {len(to_check) - amount_checked}')
        amount_checked += 1

        function = Function(phi_x, phi_y)
        checker = Checker(function)
        tests = 20
        result = checker.test(tests)
        try:
            if result[0] == CheckResult.NOT_HOLDS:
                write_not_holds(CheckInfo(function, result[0], result[1], tests), str(checker.calculate_for_point(result[1]).evalf()))
            else: 
                write_result(CheckInfo(function, result[0], result[1], tests))
        except Exception as e:
            logging.error(e.__str__())
