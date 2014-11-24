#!/usr/bin/gnuplot

filename = system("echo ${chart_filename}.png")

#PNG
set term png
set output filename
set xlabel "Reported quality"
set ylabel "Empirical quality"
set xrange [0:40]
set yrange [0:40]
stats "chart.tmp" using 1:($2 != 0 ? $2 : 1/0) name "A" 
stats "chart2.tmp" using 1:($2 != 0 ? $2 : 1/0) name "B" 
set parametric
set trange [0:40]
fx(t) = t
fy(t) = t
fAy(t) = A_slope * t + A_intercept
fBy(t) = B_slope * t + B_intercept
plot 	"chart.tmp"	using 1:($2 != 0 ? $2 : 1/0) title "Original", \
 	"chart2.tmp" 	using 1:($2 != 0 ? $2 : 1/0) title "Recalibrated", \
	fx(t),fy(t) 	lt 3 lw 2 lc rgb "blue" title "Ideal", \
 	t,fAy(t) 	lt 3 lw 1 lc rgb "red" title sprintf("Orig r=%1.4f", A_correlation), \
 	t,fBy(t) 	lt 3 lw 1 lc rgb "green" title sprintf("Recal r=%1.4f", B_correlation)

#WINDOW
#set term x11 0
#set xlabel "Reported quality"
#set ylabel "Empirical quality"
#set xrange [0:40]
#set yrange [0:40]
#set parametric
#set trange [0:40]
#fx(t) = t
#fy(t) = t
#plot 	fx(t),fy(t) 	lt 1 title "Lineal", \
#	"chart.tmp"	using 1:2 title "Original", \
# 	"chart2.tmp" 	using 1:2 title "Recalibrated" 

#pause -1
