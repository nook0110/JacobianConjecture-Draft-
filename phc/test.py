from sympy import plot_implicit, symbols, Eq
from phcpy.dimension import get_core_count
from phcpy.solver import solve
from phcpy.solutions import filter_real
from datetime import datetime

H = [ 'x**108 + 1.1*y**54 - 1.1*y;',
'y**108 + 1.1*x**54 - 1.1*x;' ]
nbcores = get_core_count()
print('Solving on', nbcores, 'cores ...')
wstart = datetime.now()
sols = solve(H, tasks=nbcores)
wstop = datetime.now()
print('  Number of solutions :', len(sols))
print('start time :', wstart)
print(' stop time :', wstop)
print('   elapsed :', wstop - wstart)