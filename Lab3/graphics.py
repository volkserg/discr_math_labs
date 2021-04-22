from solver import Solver
import PyGnuplot as gp


class Graphics:

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

    def draw_Q1_Qx0(self, a, b, g, da=0.001):
        s = Solver(a, b, g)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(s.calc_Q1())
            y2.append(s.calc_Qx0())
            a += da
            s.reload(a, b, g)
        points = []
        points.append(((x1, y1), 'tmp.dat', 'Q1(a)'))
        points.append(((x2, y2), 'tmp2.dat', 'Q*0(a)'))
        self._draw(points=points,
                   title="График зависимости Q1(a) и Q*0(a)",
                   )

    def draw_Q1_Qx1(self, a, b, g, da=0.001):
        s = Solver(a, b, g)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(s.calc_Q1())
            y2.append(s.calc_Qx1())
            a += da
            s.reload(a, b, g)
        points = []
        points.append(((x1, y1), 'tmp.dat', 'Q1(a)'))
        points.append(((x2, y2), 'tmp2.dat', 'Q*1(a)'))
        self._draw(points=points,
                   title="График зависимости Q1(a) и Q*1(a)",
                   )

    def draw_Q1_Qx0Qx1(self, a, b, g, da=0.001):
        s = Solver(a, b, g)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < b:
            x1.append(a)
            x2.append(a)
            y1.append(s.calc_Q1())
            y2.append(s.calc_Qx0() + s.calc_Qx1())
            a += da
            s.reload(a, b, g)
        points = []
        points.append(((x1, y1), 'tmp.dat', 'Q1(a)'))
        points.append(((x2, y2), 'tmp2.dat', 'Q*0(a) + Q*1(a)'))
        self._draw(points=points,
                   title="График зависимости Q1(a) и Q*0(a) + Q*1(a)",
                   )
