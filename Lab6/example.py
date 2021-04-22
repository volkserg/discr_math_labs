from solver import *

def newb(value):
    res = [0]
    al=2.02
    for i in range(1,value):
        res.append(1/i**al-1/(i+1)**al)
    return res


s = Solver(0.3, 100, newb)

# print(s.calc_pk(0))

print(s.calc_v1())
print(s.calc_v2())

print(s.calc_ing1())
print(s.calc_ing2())


print(s.calc_vhatm1(5))
print(s.calc_vhatm2(5))