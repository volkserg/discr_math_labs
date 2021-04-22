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

    def draw_N1(self, a, j, da=0.001):
        a+=da
        s1 = Solver(a, j)
        s2 = Solver(a, j+1)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < s1.inB**(-1):
            x1.append(a)
            y1.append(s1.calc_N1())
            a += da
            s1.reload(a, j)
        a = 0
        while a < s2.inB**(-1):
            x2.append(a)
            y2.append(s2.calc_N1())
            a += da
            s2.reload(a, j+1)




        points = []
        points.append(((x1, y1), 'tmp.dat', 'N1(a)'))
        points.append(((x2, y2), 'tmp2.dat', 'N1(2)(a)'))
        self._draw(points=points,
                   title="График зависимости N1(a) и N1(2)(a)",
                   )