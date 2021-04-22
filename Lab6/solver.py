import math


class Derivative:
    def __init__(self, f, h=0.0001):
        self.f = f
        self.h = float(h)

    def __call__(self, x):
        f, h = self.f, self.h
        return (f(x+h) - f(x-h))/(2*h)

class Derivative12:
    def __init__(self, f, h=0.0001):
        self.f = f
        self.h = float(h)

    def __call__(self, x, m):
        f, h = self.f, self.h
        return (f(x+h,m) - f(x-h,m))/(2*h)


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
        self.pk = {}
        self.pik = {}
        self.gk = {}
        self.vk = {}
        assert self.calc_beta(self.inA)>self.inA/(1+self.inA)

    def calc_inB(self):
        sum = 0
        for i in range(self.inf):
            sum+=i*self.b[i]
        return sum

    def calc_Bi(self, i):
        sum = 0
        for j in range(i, self.inf):
            sum += self.b[j]
        return sum

    def calc_beta(self, z):
        sum = 0
        for i in range(self.inf):
            sum += z**i*self.b[i]
        return sum

    def calc_pk(self, k):
        if self.pk.get(k) or self.pk.get(k)==0:
            return self.pk[k]
        else:
            if k==0:
                self.pk[k] = 1 - (self.inA*(1-self.calc_beta(self.inA))/self.calc_beta(self.inA))
            elif k==1:
                self.pk[k] = (1/self.calc_beta(self.inA)*(1-self.calc_beta(self.inA)))*self.calc_pk(0)
            else:
                self.pk[k] = (1/self.calc_beta(self.inA))*(1-self.calc_beta(self.inA))*(self.calc_pk(k-1)-self.calc_pik(1, k-1))
            return self.pk[k]

    def calc_pik(self, i, k):
        assert k>=1
        if self.pik.get((i, k)) or self.pik.get((i, k))==0:
            return self.pik[(i, k)]
        else:
            if i==1 and k==1:
                self.pik[(i, k)] = self.a/self.inA*self.calc_pk(0)
            elif i==1:
                self.pik[(i, k)] = self.a/self.inA*(self.calc_pk(k-1) - self.calc_pik(1, k-1))
            elif k==1:
                sum = 0
                for j in range(i, self.inf):
                    sum+=self.inA**(j-i)*self.b[j]
                self.pik[(i, k)] = self.a*(self.calc_pk(0)+self.calc_pk(1))*sum
            else:
                sum = 0
                for j in range(i, self.inf):
                    sum += self.inA**(j-i)*self.b[j]
                self.pik[(i, k)] = self.a*(self.calc_pk(k-1)-self.calc_pik(1, k-1))*sum+self.a*self.calc_pk(k)*sum
            return self.pik[(i, k)]

    def calc_gk(self, k):
        if self.gk.get(k) or self.gk.get(k)==0:
            return self.gk[k]
        else:
            if k==0:
                self.gk[k] = 0
            else:
                sum1, sum2 = 0, 0
                for i in range(1, k+1):
                    sum1 += self.inA**(i-1)*self.a*self.b[i]*self.calc_gk(k-i)
                    sum3 = 0
                    for j in range(k-i+1):
                        sum3 += self.calc_gk(j)*self.calc_gk(k-i-j)
                    sum2 += self.inA**(i-1)*self.a*self.calc_Bi(i+1)*sum3
                self.gk[k] = self.inA**k*self.b[k] + sum1 + sum2
            return self.gk[k]

    def calc_vk(self, k):
        if self.vk.get(k) or self.vk.get(k)==0:
            return self.vk[k]
        else:
            if k==0:
                self.vk[k] = 0
            else:
                sum1 = 0
                for i in range(1,k+1):
                    sum2 = 0
                    for j in range(k-i+1):
                        sum2 += self.calc_gk(j)*self.calc_vk(k-i-j)
                    sum1 += self.inA**(i-1)*self.a*self.calc_Bi(i+1)*sum2
                self.vk[k] = self.inA**(k-1)*self.b[k]+sum1
            return self.vk[k]

    def calc_vkm(self, k, m):
        if k==0:
            return 0
        elif k >= 1 and k <= m-1:
            sum1 = 0
            for i in range(1,k+1):
                sum2 = 0
                for j in range(0, k-i+1):
                    sum2 += self.calc_gk(j)*self.calc_vk(k-i-j)
                sum1 += self.inA**(i-1)*self.a*sum2
            return sum1
        else:
            sum1 = 0
            for i in range(1, m):
                sum2 = 0
                for j in range(0,k-i+1):
                    sum2 += self.calc_gk(j)*self.calc_vk(k-i-j)
                sum1 += self.inA**(i-1)*self.a*sum2
            return self.inA**(m-1) + sum1

    def calc_N(self):
        sum1 = 0
        for i in range(1, self.inf):
            sum2 = 0
            for k in range(1, self.inf):
                sum2 += k*self.calc_pik(i, k)
            sum1 += sum2
        return sum1

    def calc_gamma(self, z):
        sqr = math.sqrt((1-self.inA*z)**2*(self.inA-self.a*self.calc_beta(self.inA*z))**2 \
                        - 4*self.a*(self.inA*z-self.calc_beta(self.inA*z))*self.inA*(1-self.inA*z)*self.calc_beta(self.inA*z))
        return ((1-self.inA*z)*(self.inA-self.a*self.calc_beta(self.inA*z)) - sqr)/(2*self.a*(self.inA*z-self.calc_beta(self.inA*z)))

    def calc_ing1(self):
        df = Derivative(self.calc_gamma)
        return df(1)

    def calc_ing2(self):
        summ=0
        for k in range(self.inf):
            summ += k*self.calc_gk(k)
        return summ

    def calc_vz(self, z):
        return ((1-self.inA*z)*self.calc_beta(self.inA*z))/\
               (self.inA*(1-self.inA*z) - self.a*(self.inA*z - self.calc_beta(self.inA*z))*self.calc_gamma(z))

    def calc_v1(self):
        df = Derivative(self.calc_vz)
        return df(1)

    def calc_v2(self):
        summ = 0
        for i in range(self.inf):
            summ += i* self.calc_vk(i)
        return summ

    def calc_vzm(self, z, m):
        return ((self.inA*z)**m/self.inA)+(self.a*z/self.inA) *\
                                          (1-(self.inA*z)**(m-1))*self.calc_gamma(z)*self.calc_vz(z) /\
                                          (1-self.inA*z)

    def calc_vhatm1(self, m):
        df = Derivative12(self.calc_vzm)
        return df(1, m)

    def calc_vhatm2(self, m):
        summ = 0
        for i in range(self.inf):
            summ += i*self.calc_vkm(i, m)
        return summ