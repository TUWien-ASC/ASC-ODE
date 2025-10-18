import sys
# sys.path.append('/Users/joachim/texjs/lva/ws2324/ScientificComputing/ASC-ODE/build/mass_spring')

sys.path.append('../build/mass_spring')
sys.path.append('../build/nanoblas')

from nanoblas_impl import *
from mass_spring import *
import ngsolve.bla


mss = MassSpringSystem3d()
mss.gravity = (0,0,-9.81)

mA = mss.add (Mass(1, (1,0,0)))
mB = mss.add (Mass(2, (2,0,0)))
f1 = mss.add (Fix( (0,0,0)) )
mss.add (Spring(1, 10, (f1, mA)))
mss.add (Spring(1, 20, (mA, mB)))


print ("state = ", mss.getState())

simulate (mss, 0.1, 10)

print ("state = ", mss.getState())


simulate (mss, 0.1, 10)

print ("state = ", mss.getState())

for m in mss.masses:
    print (m.mass, m.pos)

mss.masses[0].mass = 5

for m in mss.masses:
    print (m.mass, m.pos)
