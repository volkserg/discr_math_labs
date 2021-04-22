import math


class Derivative:
    def __init__(self, f, h=0.0001):
        self.f = f
        self.h = float(h)

    def __call__(self, x):
        f, h = self.f, self.h
        return (f(x+h) - f(x-h))/(2*h)

class Derivative2:
    def __init__(self, f, h=0.0001):
        self.f = f
        self.h = float(h)

    def __call__(self, x):
        f, h = self.f, self.h
        return (f(x+h) - 2 * f(x) + f(x-h))/h**2


def C(n, k):
    return math.factorial(n)/(math.factorial(n-k)*math.factorial(k))

def b(value):
    res = [0]
    summ = 0
    for i in range(1, value):
        summ += i
    for i in range(1, value):
        res.append(i/summ)
    return res


class Solver:

    def __init__(self, a, infin=3, fun = b):
        assert 0 < a < 1
        self.a = a
        self.inf = infin
        self.b = fun(self.inf)
        self.inA = 1 - a
        self.inB = self.calc_inB()
        self.q = []
        self.p = {}
        self.ro = self.a*self.inB

    def reload(self, a, infin=3):
        assert 0 < a < 1
        self.a = a
        self.inf = infin
        self.b = b(self.inf)
        self.inA = 1 - a
        self.inB = self.calc_inB()
        self.ro = self.a*self.inB

    def calc_inB(self):
        summ = 0
        for i in range(self.inf):
            summ += i*self.b[i]
        return summ

    def calc_Bi(self, i):
        summ = 0
        for j in range(i, self.inf):
            summ += self.b[j]
        return summ

    def calc_pmn(self, m, n):
        if self.p.get((m,n)):
            return self.p.get((m,n))
        if m==0:
            if n==0:
                self.p[(m, n)] = 1 - self.a*self.calc_Bi(1)
                return self.p[(m, n)]
            elif n==1:
                self.p[(m, n)] = self.a*self.calc_Bi(1)
                return self.p[(m, n)]
            else:
                self.p[(m, n)] = 0
                return self.p[(m, n)]
        elif n==0:
            self.p[(m, n)] = ((1-self.a*self.calc_Bi(m+1))*self.calc_pmn(m-1, 0))
            return self.p[(m, n)]
        else:
            self.p[(m, n)] = ((1-self.a*self.calc_Bi(m+1))*self.calc_pmn(m-1,n) + self.a*self.calc_Bi(m+1) *\
                    self.calc_pmn(m-1, n-1))
            return self.p[(m, n)]

    def calc_Pz(self, z):
        mult = 1
        for n in range(self.inf):
            mult *= (1-self.a*(1-z)*self.calc_Bi(n+1))
        return mult

    def calc_N1(self):
        df = Derivative(self.calc_Pz)
        return df(1)

    def calc_N2(self):
        summ = 0
        for n in range(self.inf):
            summ += self.calc_Bi(n+1)
        return self.a*summ

    def calc_N3(self):
        summ = 0
        for n in range(1, self.inf):
            summ += n*self.calc_pmn(10, n)
        return summ

    def calc_D(self):
        df = Derivative(self.calc_Pz)
        df2 = Derivative2(self.calc_Pz)
        return df2(1)+df(1)-df(1)**2
