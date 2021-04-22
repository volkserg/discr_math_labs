from solver import *
from graph import *
s = Solver(0.1, 4)
# print(s.calc_wn(1)+s.calc_wn(0))
# print((1-s.ro)/s.inA)
# # print(s.b)
#
#
# print(s.calc_Pwz(2))
# print(s.calc_betaz(s.a*2+s.inA))

# print(s.calc_N1())
# print(s.calc_N2())
# print(s.calc_N3())
# print(s.calc_N4())
#
#
# print(s.calc_Nx1())
# print(s.calc_Nx2())
# print(s.calc_Nx3())



# print(s.calc_M1())
# print(s.calc_M2())
# print(s.calc_M3())

# print(s.check_Pxz(0.9))

print(s.calc_pplusi(0)+s.calc_pplusi(1))
print((1-s.ro)/s.calc_pwn(0))


# g = Graphics()
# g.draw_N1_Nx1(0, 4)