import math


class Derivative:
    def __init__(self, f, h=0.0001):
        self.f = f
        self.h = float(h)

    def __call__(self, x):
        f, h = self.f, self.h
        return (f(x+h) - f(x-h))/(2*h)


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

    def __init__(self, a, infin=3):
        assert 0 < a < 1
        self.a = a
        self.inf = infin
        self.b = b(self.inf)
        self.inA = 1 - a
        self.inB = self.calc_inB()
        self.q = []
        self.ro = self.a*self.inB

    def reload(self, a, infin=3):
        assert 0 < a < 1
        self.a = a
        self.inf = infin
        self.b = b(self.inf)
        self.inA = 1 - a
        self.inB = self.calc_inB()
        self.q = []
        self.ro = self.a*self.inB

    def calc_inB(self):
        res = 0
        for i in range(self.inf):
            res += i * self.b[i]
        return res

    def calc_Bi(self, i):
        res = 0
        for j in range(i,self.inf):
            res += self.b[j]
        return res

    def calc_pwn(self, n):
        assert n >= 0
        res = 0
        for t in range(n, self.inf):
            res += C(t, n)*self.a**n*self.inA**(t-n)*self.b[t]
        return res

    def calc_Pwz(self, z):
        res = 0
        for n in range(self.inf):
            res += z**n*self.calc_pwn(n)
        return res

    def calc_Pwi(self, i):
        res = 0
        for j in range(i, self.inf):
            res += self.calc_pwn(j)
        return res

    def calc_betaz(self, z):
        summ = 0
        for n in range(self.inf):
            summ+= z**n*self.b[n]
        return summ

    def calc_qplusi(self, i):
        while len(self.q)<=i:
            self.q.append(None)
        if self.q[i]:
            return self.q[i]
        else:
            tmp = 0
            for j in range(1, i):
                tmp += self.calc_qplusi(j)*self.calc_Pwi(i-j+1)
            self.q[i] = (self.calc_Pwi(i)+tmp)/self.calc_pwn(0)
            return self.q[i]

    def calc_pplusi(self, i):
        if i == 0:
            summ = 0
            for j in range(1,self.inf):
                summ+=self.calc_qplusi(j)
            return 1/(1+summ)
        else:
            return self.calc_pplusi(0)*self.calc_qplusi(i)

    def calc_Pplusz(self, z):
        # assert 0 < z < 1
        return (1-z)*self.calc_Pwz(z)*(1-self.a*self.inB)/(self.calc_Pwz(z)-z)

    def calc_phat(self, i):
        assert i>=1
        if i == 1:
            return self.calc_pplusi(0) + self.calc_pplusi(1)
        else:
            return self.calc_pplusi(i)

    def calc_pi(self, i):
        if i == 0:
            return self.calc_pplusi(0)
        else:
            summ = 0
            for n in range(0, self.inf):
                for j in range(max(1,i-n), i+1):
                    summ += self.calc_phat(j)*self.calc_Bi(n+1)*C(n,i-j)*self.a**(i-j)*self.inA**(n-i+j)
            return self.a*summ

    def check_pi(self, i):
        summ1, summ2 = 0, 0
        for n in range(i-1, self.inf):
            summ1 += self.calc_Bi(n+1)*C(n,i-1)*self.a**(i-1)*self.inA**(n-i+1)
        for j in range(1,i+1):
            summ3 = 0
            for n in range(i-j,self.inf):
                self.calc_Bi(n+1)*C(n,i-j)*self.a**(i-j)*self.inA**(n-i+j)
            summ2+=self.calc_pplusi(i)*summ3
        return self.a*(self.calc_pplusi(0)*summ1+summ2)

    def calc_Pz(self, z):
        return (1 - self.calc_Pwz(z))*((self.calc_Pplusz(z)-self.calc_pplusi(0)) \
                                       +z*self.calc_pplusi(0))/(1-z)+self.calc_pplusi(0)

    def calc_pxi(self, i):
        if i == 0:
            return self.calc_pplusi(0)/self.inA
        else:
            summ1, summ2, summ3, summ4 = 0, 0, 0, 0
            for n in range(i,self.inf):
                summ1 += self.calc_Bi(n+1)*C(n-1,i-1)*self.a**(i-1)*self.inA**(n-i)

            for j in range(1,i+1):
                summ = 0
                for n in range(i-j+1, self.inf):
                    summ+= self.calc_Bi(n+1)*C(n-1,i-j)*self.a**(i-j)*self.inA**(n-i+j-1)
                summ2+=self.calc_pplusi(j)*summ

            for n in range(i+1, self.inf):
                summ3 += self.b[n]*C(n-1,i)*self.a**i*self.inA**(n-i-1)

            for j in range(1,i+2):
                summ = 0
                for n in range(i-j+2, self.inf):
                    summ += self.b[n]*C(n-1,i-j+1)*self.a**(i-j+1)*self.inA**(n-i+j-2)
                summ4 += self.calc_pplusi(j)*summ

            return self.a*(self.calc_pplusi(0)*summ1+summ2)+self.a*self.calc_pplusi(0)*summ3+self.a*summ4

    def calc_Pxz(self, z):
        res = 0
        for i in range(self.inf):
            res += z**i*self.calc_pxi(i)
        return res

    def check_Pxz(self, z):
        return self.calc_pplusi(0)/self.inA + \
    self.a*self.calc_pplusi(0)*z*\
               ((self.a*z+self.inA)-self.calc_betaz(self.a*z+self.inA))/\
               ((self.a*z+self.inA)*(1-(self.a*z+self.inA)))+\
    self.a*(self.calc_Pplusz(z)-self.calc_pplusi(0))*\
               ((self.a*z+self.inA)-self.calc_betaz(self.a*z+self.inA))/\
               ((self.a*z+self.inA)*(1-(self.a*z+self.inA)))+\
    self.a*self.calc_pplusi(0)*((self.calc_betaz(self.a*z+self.inA)/(self.a*z+self.inA))-\
                                self.calc_betaz(self.inA)/self.inA)+\
    self.a/z*(self.calc_Pplusz(z)-self.calc_pplusi(0))*\
               (self.calc_betaz(self.a*z+self.inA)/(self.a*z+self.inA)-self.calc_betaz(self.inA)/self.inA)+\
    self.a/self.inA*(self.calc_Pplusz(z)-self.calc_pplusi(0)-self.calc_pplusi(1)*z)*\
               (self.calc_betaz(self.inA)-self.b[0])

    def calc_wn(self, n, w=[]):
        if n == 0:
            return 1 - self.ro
        elif n==1:
            return self.a*self.calc_wn(0)/self.inA
        else:
            while len(w)<=n:
                w.append(None)
            if w[n]:
                return w[n]
            else:
                summ = 0
                for i in range(n-2):
                    summ+=self.calc_wn(i+1, w)*self.a\
                          *self.b[n-1-i]
                w[n] = 1/self.inA*(self.calc_wn(n-1, w)*(1-self.a*self.b[1]) \
                                   - self.calc_wn(0)*self.a*self.b[n-1] - summ)
                return w[n]

    def calc_whatn(self, n):
        if n ==0:
            return self.calc_wn(0)+self.calc_wn(1)
        else:
            return self.calc_wn(n+1)

    def calc_wz(self, z):
        return (z-1)*(self.inA+self.a*self.calc_betaz(z))*self.calc_wn(0)/(z-self.inA-self.a*self.calc_betaz(z))

    def calc_N1(self):
        df = Derivative(self.calc_Pplusz)
        return df(1)

    def calc_N2(self):
        summ = 0
        for i in range(self.inf):
            summ += i*self.calc_pi(i)
        return summ

    def calc_N3(self):
        df = Derivative(self.calc_Pz)
        return df(1)

    def calc_N4(self):
        summ = 0
        for i in range(self.inf):
            summ+=i*self.calc_pi(i)
        return summ

    def calc_Nx1(self):
        df = Derivative(self.calc_Pxz)
        return df(1)

    def calc_Nx2(self):
        df = Derivative(self.check_Pxz)
        return df(1)

    def calc_Nx3(self):
        summ = 0
        for i in range(self.inf):
            summ += i* self.calc_pxi(i)
        return summ

    def calc_M1(self):
        summ = 0
        for i in range(self.inf):
            summ += i*self.calc_wn(i)
        return summ

    def calc_M2(self):
        df = Derivative(self.calc_wz)
        return df(1)

    def calc_M3(self):
        summ = 0
        for i in range(self.inf):
            summ += i* self.calc_whatn(i)
        return summ