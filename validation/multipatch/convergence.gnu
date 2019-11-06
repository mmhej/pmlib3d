reset

set terminal postscript eps enhanced

set output 'convergence.eps'
set size 0.65,0.65
set pointsize 1
## Format the axes
set log x
set log y
set format x "10^{%T}"
set format y "10^{%T}"
## Set labels - use X & Y and psfrag
set yl 'rms error'
set xl 'h'
## Set plotting range. tran is for parametric plots (used for vertical lines)
set xran [ 1e-4 : 0.1 ]
set yran [ 1e-14 : 1e0 ]
## Set tics
#set xtics 0.2 
#set ytics 2

## Plot command (visual/X11)
plot  \
  1e2*x**2 with line linetype 2 lc rgb 'black' lw 1 notitle,\
  1e4*x**4 with line linetype 2 lc rgb 'black' lw 1 notitle,\
  1e6*x**6 with line linetype 2 lc rgb 'black' lw 1 notitle,\
  1e8*x**8 with line linetype 2 lc rgb 'black' lw 1 notitle,\
  1e10*x**10 with line linetype 2 lc rgb 'black' lw 1 notitle,\
  'convergence.dat' u ($1/0.5):2 with line linetype 1 lc rgb 'blue' lw 1 notitle



