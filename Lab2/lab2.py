import math
import PyGnuplot as gp


def C(n, k):
    return math.factorial(n)/(math.factorial(n-k)*math.factorial(k))


class Solver:

    def __init__(self, a, b, n):
        assert a > 0 and a < 1 and b > 0 and b < 1 and n >= 1
        self.a = a
        self.b = b
        self.inA = 1 - a
        self.inB = 1 - b
        self.n = n
        self.count = 0
        self.q = [None for i in range(n+1)]
        self.q[n] = 1

    def reload(self, a, b, n):
        assert a > 0 and a < 1 and b > 0 and b < 1 and n >= 1
        self.a = a
        self.b = b
        self.inA = 1 - a
        self.inB = 1 - b
        self.n = n
        self.q = [None for i in range(n+1)]
        self.q[n] = 1
    # -------------------------1st----------------------------
    def calc_pij(self, i, j):
        if j == 0:
            return self.inA*self.b**i
        elif j == i+1 and i < self.n:
            return self.a*self.inB**i
        elif i == self.n and j == self.n:
            return self.inB**self.n + self.a*C(self.n, 1)*self.inB**(self.n-1)*self.b
        else:
            return self.inA*C(i, j)*(self.inB**j)*self.b**(i-j) + self.a * C(i, j-1) * (self.inB**(j-1))*self.b**(i-j+1)

    def calc_Pij(self, i, j):
        res = 0
        assert i >= 1 and i <= self.n and j >= 0 and j <= i-1
        for k in range(0, j+1):
            res += self.calc_pij(i, k)
        return res

    def calc_qi(self, i):
        if self.q[i]:
            return self.q[i]
        else:
            tmp = 0
            for j in range(self.n, i, -1):
                tmp += self.calc_qi(j)*self.calc_Pij(j, i)
            self.q[i] = tmp/self.calc_pij(i, i+1)
            return self.q[i]

    def calc_pi(self, i):
        assert i >= 0 and i <= self.n
        if i == self.n:
            tmp = 0
            for k in range(0, self.n):
                tmp += self.calc_qi(k)
            return 1/(1+tmp)
        else:
            return self.calc_pi(self.n) * self.calc_qi(i)

    def calc_pxji(self, j, i):
        return C(j, i) * self.inB**i * self.b**(j-i)

    def calc_pxi(self, i):
        sum = 0
        for j in range(i, self.n+1):
            sum += self.calc_pi(j)*self.calc_pxji(j, i)
        return sum

    # -------------------------2nd----------------------------

    def calc_special_pi(self):
        return self.calc_pxi(self.n)

    def calc_N(self):
        res = 0
        for i in range(self.n+1):
            res += i * self.calc_pi(i)
        return res

    def calc_Nx(self):
        res = 0
        for i in range(self.n+1):
            res += i * self.calc_pxi(i)
        return res

    def calc_VarN(self):
        res = 0
        for i in range(self.n+1):
            res += (i - self.calc_N())**2*self.calc_pi(i)
        return res

    def calc_VarNx(self):
        res = 0
        for i in range(self.n+1):
            res += (i - self.calc_Nx())**2*self.calc_pxi(i)
        return res


class Graphics():

    def _draw(self,
             points, # ([x1,y1],filename,functionname), ...
             xl='Значения a',
             yl='Значения функции',
             title='заголовок',
             yrange='[0:5]',
             xrange='[-1:1]',
             out_file='file.pdf'):
        gp.c('set xlabel "' + xl + '"')
        gp.c('set ylabel "' + yl + '"')
        gp.c('set title "' + title + '"')
        gp.c('set yrange ' + yrange)
        gp.c('set xrange ' + xrange)
        plotstr = 'plot '
        for q in points:
            gp.s([q[0][0], q[0][1]], filename=q[1])
            plotstr += '"' + q[1] + '" u 1:2 w l title "' + q[2] + '", '
        plotstr = plotstr.strip(', ')
        gp.c(plotstr)
        # print(plotstr)
        # gp.pdf("out.pdf")

    def draw_N_Nx(self, a, b, n, da=0.001):
        s = Solver(a, b, n)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < 1-da:
            x1.append(a)
            x2.append(a)
            y1.append(s.calc_N())
            y2.append(s.calc_Nx())
            a += da
            s.reload(a, b, n)
        points = []
        points.append(((x1, y1), 'tmp.dat', 'N(a)'))
        points.append(((x2, y2), 'tmp2.dat', 'N*(a)'))
        self._draw(points=points,
                   title="График зависимости N(a) и N*(a)",
                   )

    def draw_VarN_VarNx(self, a, b, n, da=0.001):
        s = Solver(a, b, n)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < 1-da:
            print(a)
            x1.append(a)
            x2.append(a)
            y1.append(s.calc_VarN())
            y2.append(s.calc_VarNx())
            a += da
            s.reload(a, b, n)
        points = []
        points.append(((x1, y1), 'tmp.dat', 'VarN(a)'))
        points.append(((x2, y2), 'tmp2.dat', 'VarN*(a)'))
        self._draw(points=points,
                   title="График зависимости VarN(a) и VarN*(a)",
                   )

    def draw_special_pi(self, a, b, n, da=0.001):
        s = Solver(a, b, n)
        x1 = []
        y1 = []
        while a < 1-da:
            print(a)
            x1.append(a)
            y1.append(s.calc_special_pi())
            a += da
            s.reload(a, b, n)
        points = []
        points.append(((x1, y1), 'tmp.dat', 'p*n(a)'))
        self._draw(points=points,
                   title="График зависимости VarN(a) и VarN*(a)",
                   )
