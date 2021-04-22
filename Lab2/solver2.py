from lab2 import C

class Solver2:

    def __init__(self, a, b, n, r):
        assert 0 < a < 1 and 0 < b < 1 <= n and r >= 1
        self.a = a
        self.b = b
        self.inA = 1 - a
        self.inB = 1 - b
        self.n = n
        self.r = r
        self.R = n + r
        self.count = 0
        self.q = [None for i in range(n+r+1)]
        self.q[n+r] = 1

    # -------------------------1st----------------------------

    def calc_pij(self, i, j):

        if i == j == self.R:
            return self.inB**self.n + self.a*self.n*self.inB**(self.n-1)*self.b
        elif j == i+1 and self.n <= i <= self.R-1:
            return self.a*self.inB**self.n
        elif i - self.n + 1 <= j <= min(i, self.R-1) and self.n <= i <= self.R:
            return self.inA*C(self.n, i-j)*self.inB**(self.n-i+j)*self.b**(i-j) + self.a*C(self.n, i-j+1)*self.inB**(self.n-i+j-1)*self.b**(i-j+1)
        elif self.n <= i <= self.R and j == i - self.n:
            return self.inA*self.b**self.n
        elif 0 <= i <= self.n - 1 and j == i + 1:
            return self.a*self.inB**i
        elif self.n - 1 >= i >= 0 == j:
            return self.inA*self.b**i
        elif 1 <= i <= self.n - 1 and 1 <= j <= i:
            return self.inA*C(i, j)*self.inB**j*self.b**(i-j) + self.a*C(i, j-1)*self.inB**(j-1)*self.b**(i-j+1)

    def calc_qi(self, i):
        if self.q[i]:
            if self.q[i-1]:
                return self.q[i]
            else:
                if self.r >= i >= 1:
                    summ = 0
                    for j in range(i+1, i+self.n+1):
                        summ += self.calc_pij(j, i)*self.calc_qi(j)
                    self.q[i-1] = ((1 - self.calc_pij(i, i))*self.q[i] - summ) / self.calc_pij(i-1, i)
                elif self.R >= i >= self.r+1:
                    summ = 0
                    for j in range(i+1, self.R+1):
                        summ += self.calc_pij(j, i)*self.calc_qi(j)
                    self.q[i-1] = ((1 - self.calc_pij(i, i))*self.q[i] - summ) / self.calc_pij(i-1, i)
                return self.q[i]
        else:
            self.calc_qi(i+1)
            return self.calc_qi(i)

    def calc_pi(self, i):
        if i == self.R:
            summ = 0
            for k in range(0, self.R):
                summ += self.calc_qi(k)
            return 1 / (1 + summ)
        else:
            return self.calc_pi(self.R)*self.calc_qi(i)

    def calc_pxij(self, i, j):
        if self.n <= i <= self.R and i - self.n <= j <= i:
            return C(self.n, i-j)*self.inB**(self.n-i+j)*self.b**(i-j)
        elif 0 <= i <= self.n - 1 and 0 <= j <= i:
            return C(i, j)*self.inB**j*self.b**(i-j)

    def calc_pxi(self, i):
        if self.r+1 <= i <= self.R:
            summ = 0
            for j in range(i, self.R+1):
                summ += self.calc_pi(j) * self.calc_pxij(j, i)
        elif 0 <= i <= self.r:
            summ = 0
            for j in range(i, i+self.n+1):
                summ += self.calc_pi(j) * self.calc_pxij(j, i)
        else:
            summ = 'errrrrrror'
        return summ

    def calc_special_pi(self):
        return self.calc_pxi(self.R)

    def calc_wxij(self, i, j):
        if i == 1:
            if 0 <= j <= self.n:
                return C(self.n, j)*self.b**j*self.inB**(self.n-j)
            else:
                return 0
        else:
            if 0 <= j <= i*self.n and i >= 2:
                return C(i*self.n, j)*self.b**j*self.inB**(self.n*i-j)
            else:
                return 0

    def calc_wij(self,i, j):
        if i == 1:
            summ = 0
            for k in range(j+1, 1000):
                summ += self.calc_wxij(1, k)
            return summ
        else:
            summ = 0
            for k in range(j+1):
                summ += self.calc_wxij(i-1, k) * self.calc_wij(1, j-k)
            return summ

    def calc_wi(self, i):
        if i == 0:
            summ = 0
            for j in range(self.n):
                summ += self.calc_pxi(j)
            return summ / (1 - self.calc_special_pi())
        else:
            summ = 0
            for j in range(self.n, self.R):
                summ += self.calc_pxi(j) * self.calc_wij(i, j-self.n)
            return summ / (1 - self.calc_special_pi())

    # -------------------------2nd----------------------------
    def calc_N(self):
        summ = 0
        for i in range(self.R+1):
            summ += i * self.calc_pi(i)
        return summ

    def calc_Nx(self):
        summ = 0
        for i in range(self.R+1):
            summ += i* self.calc_pxi(i)
        return summ

    def calc_VarN(self):
        summ = 0
        for i in range(self.R+1):
            summ += (i - self.calc_N)**2*self.calc_pi(i)
        return summ

    def calc_VarNx(self):
        summ = 0
        for i in range(self.R+1):
            summ += (i + self.calc_Nx())**2*self.calc_pi(i)
        return summ

    def calc_w(self):
        summ = 0
        for i in range(1000):
            summ += i * self.calc_wi(i)
        return summ

    def calc_Dw(self):
        summ = 0
        w = self.calc_w()
        for i in range(1000):
            summ += (i - w)**2*self.calc_wi(i)
        return summ
