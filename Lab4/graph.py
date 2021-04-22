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

    def draw_N1_Nx1(self, a, j, da=0.001):
        a+=da
        s = Solver(a, j)
        x1 = []
        y1 = []
        x2 = []
        y2 = []
        while a < s.inB**(-1):
            x1.append(a)
            x2.append(a)
            y1.append(s.calc_N1())
            y2.append(s.calc_Nx1())
            a += da
            s.reload(a, j)
        points = []
        points.append(((x1, y1), 'tmp.dat', 'N1(a)'))
        points.append(((x2, y2), 'tmp2.dat', 'N*1(a)'))
        self._draw(points=points,
                   title="График зависимости N1(a) и N*1(a)",
                   )