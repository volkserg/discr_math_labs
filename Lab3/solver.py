class Derivative:
    def __init__(self, f, h=0.0001):
        self.f = f
        self.h = float(h)

    def __call__(self, x):
        f, h = self.f, self.h
        return (f(x+h) - f(x))/h


class Solver:

    def __init__(self, a, b, g, infin=10):
        assert 0 < a < 1 and 0 < b < 1 and 0 < g < 1
        self.a = a
        self.b = b
        self.g = g
        self.inA = 1 - a
        self.inB = 1 - b
        self.inG = 1 - g
        self.q0 = [1]
        self.q1 = [a / (self.inA * b)]
        self.rho = a / b
        self.P1 = self.rho
        self.P0 = 1 - self.rho
        self.inf = infin

    def reload(self, a, b, g, infin=10):
        assert 0 < a < 1 and 0 < b < 1 and 0 < g < 1
        self.a = a
        self.b = b
        self.g = g
        self.inA = 1 - a
        self.inB = 1 - b
        self.inG = 1 - g
        self.q0 = [1]
        self.q1 = [a / (self.inA * b)]
        self.rho = a / b
        self.P1 = self.rho
        self.P0 = 1 - self.rho
        self.inf = infin
    # -------------------------1st-----------------

    def calc_q0(self, i):
        assert i >= 0
        if len(self.q0) <= i:
            # for k in range(len(self.q0), i):
            while len(self.q0) < i+ 1:
                self.q0.append(None)
        if self.q0[i]:
            return self.q0[i]
        else:
            self.q0[i] = self.a * self.inB * self.inG ** i * self.calc_q1(i - 1) / (1 - self.inG ** i)
            return self.q0[i]

    def calc_q1(self, i):
        assert i >= 0
        if len(self.q1) <= i:
            # for k in range(len(self.q1), i):
            while len(self.q1) < i + 1:
                self.q1.append(None)
        if self.q1[i]:
            return self.q1[i]
        else:
            self.q1[i] = self.a * self.inB * (1 - self.inA * self.inG ** i) * self.calc_q1(i - 1) / (
                        self.inA * self.b * (1 - self.inG ** i))
            return self.q1[i]

    def calc_p00(self):
        summ1, summ2 = 0, 0
        for i in range(self.inf):
            summ1 += self.calc_q0(i)
        for i in range(self.inf):
            summ2 += self.calc_q1(i)
        return 1 / (summ1 + summ2)

    def calc_p0i(self, i):
        return self.calc_p00() * self.calc_q0(i)

    def calc_p1i(self, i):
        return self.calc_p00() * self.calc_q1(i)

    def calc_p00_check(self):
        mult = 1
        for i in range(1, self.inf):
            mult *= self.inA * (self.b - self.a * self.inB * self.inG ** i) / (
                        self.inA * self.b - self.a * self.inB * self.inG ** i)
        return (1 - self.rho) * (1 / mult)

    def calc_pxi(self, i):
        if i == 0:
            return self.P0 + self.P1 * self.b
        else:
            return self.calc_p1i(i - 1) * self.inB

    def calc_vk0(self, k, n):
        assert n >= 0 and k >= 1
        if k == 1:
            return self.inA * (1 - self.inG ** (n + 1) / (n + 1))
        elif n == 0 and k >= 2:
            return self.a * self.calc_vk1(k - 1, 0) + self.inA * self.inG * self.calc_vk0(k - 1, 0)
        else:
            return self.a * self.calc_vk1(k - 1, n) + self.inA*(
                self.inG ** (n + 1) * self.calc_vk0(k - 1, n) + n * (1 - self.inG ** (n + 1)) * self.calc_vk1(k - 1,
                                                                                                              n - 1) / (
                            n + 1))

    def calc_vk1(self, k, n):
        assert n >= 0 and k >= 1
        if k == 1:
            return self.inA * self.b * (1 - self.inG ** (n + 1) / (n + 1))
        elif n == 0 and k >= 2:
            return self.a * (self.inB * self.calc_vk1(k - 1, 1) + self.b * self.calc_vk1(k - 1, 0)) + \
                   self.inA * (self.inB * self.calc_vk1(k - 1, 0) + self.b * self.inG * self.calc_vk0(k - 1, 0))
        else:
            return self.a * (self.inB * self.calc_vk1(k - 1, n + 1) + self.b * self.calc_vk1(k - 1, n)) + \
                   self.inA * (self.inB * self.calc_vk1(k - 1, n) + self.b * self.inG ** (n + 1) * self.calc_vk0(k - 1,
                                                                                                                 n) + \
                               n * self.b * (1 - self.inG ** (n + 1)) * self.calc_vk1(k - 1, n - 1) / (n + 1))

    # --------------------2nd-------------------------

    def calc_P0(self, z):
        # assert 0 < z <= 1
        mult = 1
        for i in range(1,self.inf):
            mult *= 1 + self.a**2*self.inB*self.inG**i*z/(self.inA*self.b - self.a*self.inB*self.inG**i*z)
        return self.calc_p00() * mult

    def calc_P1(self, z):
        return self.a * self.calc_P0(z) / (self.inA*self.b -self.a*self.inB*z)

    def calc_Q1(self):
        df1 = Derivative(self.calc_P0)
        df2 = Derivative(self.calc_P1)
        return df1(1) + df2(1)

    def calc_Q2(self):
        df = Derivative(self.calc_P0)
        return self.a**2*self.inB/(self.b*(self.b-self.a)) + self.b*df(1)/(self.b-self.a)

    def calc_Q3(self):
        sum1 = 0
        sum2 = 0
        for i in range(1,self.inf):
            sum1 += i*self.calc_p0i(i)
            sum2 += i*self.calc_p1i(i)
        return sum1 + sum2

    def calc_Qx0(self):
        df1 = Derivative(self.calc_P0)
        df2 = Derivative(self.calc_P1)
        return (df1(1) + self.b*df2(1))/self.calc_pxi(0)

    def calc_Qx1(self):
        df = Derivative(self.calc_P1)
        return self.inB*df(1)/(1-self.calc_pxi(0))

    def calc_wk(self, k):
        assert k >= 0
        if k == 0:
            return self.calc_pxi(0)
        else:
            sum1 = 0
            for i in range(1,self.inf):
                sum1 += self.calc_pxi(i)*self.calc_vk1(k, i-1)
            return sum1

    def calc_w(self):
        sum1 = 0
        for i in range(1, self.inf):
            sum2 = 0
            for k in range(1, self.inf):
                sum2 += k*self.calc_vk1(k, i-1)
            sum1 += self.calc_pxi(i) * sum2
        return sum1
