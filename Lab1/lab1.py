import PyGnuplot as gp
# import numpy as np
import matplotlib.pyplot as plt

class Derivative2:
    def __init__(self, f, h=0.0001):
        self.f = f
        self.h = float(h)

    def __call__(self, x):
        f, h = self.f, self.h
        return (f(x+h) - 2 * f(x) + f(x-h))/h**2


class Task1:

    def __init__(self, a, b):
        assert (a > 0 and a < 1 and b > 0 and b < 1)
        self.a = a
        self.b = b
        self.inA = 1 - a
        self. inB = 1 - b
        self.p = self.a/self.b
        self.delta = self.a * self.inB / (self.inA * self.b)

    def reload(self, a, b):
        assert (a > 0 and a < 1 and b > 0 and b < 1)
        self.a = a
        self.b = b
        self.inA = 1 - a
        self. inB = 1 - b
        self.p = self.a/self.b
        self.delta = self.a * self.inB / (self.inA * self.b)

    def get_a(self):
        return self.a

    def get_b(self):
        return self.b

    def increase_a(self, n):
        self.a += n

# ---------------------------------1st-----------------------------
    # 1-2
    def calculate_p(self, i):
        assert (i >= 0)
        p0 = 1 - self.p
        if i == 0:
            return p0
        elif i > 0:
            return (1/self.inB) * (self.delta ** i) * p0

    # 3-4
    def calculate_px(self, i):
        assert (i >= 0)
        if i == 0:
            return 1 - self.delta
        elif i >= 1:
            return (1 - self.delta)*(self.delta**i)

    # 5-6
    def calculate_vk(self, k):
        assert (k >= 1)
        return (1 - self.inB / self.inA) * (self.inB / self.inA)**(k - 1)

    # 7-8
    def calculate_wk(self, k):
        assert (k >= 0)
        if k == 0:
            return self.calculate_px(0)
        elif k >= 1:
            return (1 - self.calculate_px(0)) * self.calculate_vk(k)
    # ---------------------------------2nd-----------------------------

    # 9
    def calculate_n(self):
        return (self.inA * self.p) / (1 - self.p)

    # 10
    def calculate_q(self):
        return (self.inB * (self.p**2)) / (1 - self.p)

    # 11
    def calculate_nx(self):
        return (self.inB * self.p) / (1 - self.p)

    # 12
    def calculate_qx(self):
        return (self.inB * (self.p**2)) / (self.inA * (1 - self.p))

    # 13
    def calculate_P(self, z):
        return self.calculate_p(0) * (self.inA * self.b + z * self.a * self.b) / (self.inA * self.b - z * self.a * self.inB)

    def calculate_VarN(self):
        df = Derivative2(self.calculate_P)
        return df(1) + self.calculate_n() - self.calculate_n()**2

    # 14
    def calculate_Px(self, z):
        return (1 - self.delta) / (1 - self.delta * z)

    def calculate_VarNx(self):
        df = Derivative2(self.calculate_Px)
        return df(1) + self.calculate_nx() - self.calculate_nx()**2
    # 15
    def calculate_w(self):
        return (self.a * self.inB) / (self.b * (self.inA - self.inB))

    # 16
    def calculate_v(self):
        return self.calculate_w() + 1 / self.b

    # 17
    def calculate_wz(self, z):
        return self.calculate_px(0) + (1 - self.calculate_px(0)) * ((self.inA - self.inB) * z) / (self.inA - self.inB*z)

    def calculate_Dw(self):
        df = Derivative2(self.calculate_wz)
        return df(1) + self.calculate_w() - self.calculate_w()**2
    # 18
    def calculate_fi(self, z):
        return (self.inA - self.inB) * z / (self.inA - self.inB * z)

    def calculate_Dv(self):
        df = Derivative2(self.calculate_fi)
        return df(1) + self.calculate_v() - self.calculate_v()**2
# ---------------------------------3rd-----------------------------


class Graphics():

    def draw_Na(self, a, b, da):
        t = Task1(a, b)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(t.calculate_n())
            y2.append(t.calculate_nx())
            a += da
            t.reload(a, b)
        gp.s([x1, y1], filename='tmp1.dat')
        gp.s([x2, y2], filename='tmp2.dat')
        gp.c('set xlabel "Значение а"')
        gp.c('set ylabel "Значение функции"')
        gp.c('set title "Графики зависимостей показателей N и N* от значений a"')
        gp.c('set yrange [0:5]')
        gp.c('plot "tmp1.dat" u 1:2 w i title "N(a)", "tmp2.dat" u 1:2 w i title "N*(a)"')
        gp.pdf("NA.pdf")

    def draw_VarN(self, a, b, da):
        t = Task1(a, b)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(t.calculate_VarN())
            y2.append(t.calculate_VarNx())
            a += da
            t.reload(a, b)
        gp.s([x1, y1], filename='tmp1.dat')
        gp.s([x2, y2], filename='tmp2.dat')
        gp.c('set xlabel "Значение а"')
        gp.c('set ylabel "Значение функции"')
        gp.c('set title "Графики зависимостей показателей VarN и VarN* от значений a"')
        gp.c('set yrange [0:5]')
        gp.c('plot "tmp1.dat" u 1:2 w i title "VarN(a)", "tmp2.dat" u 1:2 w i title "VarN*(a)"')
        gp.pdf("VarNa.pdf")

    def draw_Q(self, a, b, da):
        t = Task1(a, b)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(t.calculate_q())
            y2.append(t.calculate_qx())
            a += da
            t.reload(a, b)
        gp.s([x1, y1], filename='tmp1.dat')
        gp.s([x2, y2], filename='tmp2.dat')
        gp.c('set xlabel "Значение а"')
        gp.c('set ylabel "Значение функции"')
        gp.c('set title "Графики зависимостей показателей Q и Q* от значений a"')
        gp.c('set yrange [0:5]')
        gp.c('plot "tmp1.dat" u 1:2 w i title "VarN(a)", "tmp2.dat" u 1:2 w i title "VarN*(a)"')
        gp.pdf("Q.pdf")

    def draw_vw(self, a, b, da):
        t = Task1(a, b)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(t.calculate_w())
            y2.append(t.calculate_v())
            a += da
            t.reload(a, b)
        gp.s([x1, y1], filename='tmp1.dat')
        gp.s([x2, y2], filename='tmp2.dat')
        gp.c('set xlabel "Значение а"')
        gp.c('set ylabel "Значение функции"')
        gp.c('set title "Графики зависимостей показателей w и v от значений a"')
        gp.c('set yrange [0:5]')
        gp.c('plot "tmp1.dat" u 1:2 w l title "w(a)", "tmp2.dat" u 1:2 w l title "v(a)"')
        gp.pdf("vw.pdf")

    def draw_D(self, a, b, da):
        t = Task1(a, b)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(t.calculate_Dw())
            y2.append(t.calculate_Dv())
            a += da
            t.reload(a, b)
        gp.s([x1, y1], filename='tmp1.dat')
        gp.s([x2, y2], filename='tmp2.dat')
        gp.c('set xlabel "Значение а"')
        gp.c('set ylabel "Значение функции"')
        gp.c('set title "Графики зависимостей показателей Dw и Dv от значений a"')
        gp.c('set yrange [0:5]')
        gp.c('plot "tmp1.dat" u 1:2 w l title "Dw(a)", "tmp2.dat" u 1:2 w l title "Dv(a)"')
        gp.pdf("D.pdf")
# g = Graphics()
# g.draw_D(0.1, 0.7, 0.001)
# t = Task1(0.3,0.4)
# print(t.calculate_VarN())