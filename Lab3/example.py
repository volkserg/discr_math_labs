from solver import Solver
from graphics import Graphics
s = Solver(0.5, 0.7, 0.9)


print(s.calc_w()*s.a)
print(s.calc_Q1())
print(s.calc_Q2())
print(s.calc_Q3())


# g = Graphics()
# g.draw_Q1_Qx0Qx1(0.01, 0.9, 0.5)