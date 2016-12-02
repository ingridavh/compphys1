#Computes derivatives of wavefunctions
#to be used in c++ program

from sympy import symbols, diff, exp, sqrt, factor, Symbol, printing
x, y, z, Z = symbols('x y z Z')
r = sqrt(x*x+y*y+z*z)
phi= (Z*r-2)*exp(-Z*r/2)
#Creates a symbolic equivalent of r
R = Symbol('r')
#print latex and c++ code
print printing.latex(diff(phi, x).factor().subs(r,R))
print printing.ccode(diff(phi, x).factor().subs(r,R))
#second derivative
print printing.latex((diff(diff(phi,x), x) + diff(diff(phi, x), x) + diff(diff(phi, z), z)).factor().collect(Z).subs(r,R).subs(r**2, R**2).factor())
