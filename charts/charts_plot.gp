#!/usr/bin/gnuplot

#PNG
set term png
set output "total.png"
set xlabel "Threads"
set ylabel "Time (s)"
set title "Time to process"
set yrange [0:] writeback
plot 	"<awk '{if($1 == 0){print $2 \" \" $3}}' cook.dat" w lp t "Total" 
set yrange [*:*]

set output "framework.png"
set xlabel "Threads"
set ylabel "Time (s)"
set title "Framework pipeline"
plot 	"<awk '{if($1 == 2){print $2 \" \" $3}}' cook.dat" w lp t "Read", \
 	"<awk '{if($1 == 3){print $2 \" \" $3}}' cook.dat" w lp t "Process", \
 	"<awk '{if($1 == 4){print $2 \" \" $3}}' cook.dat" w lp t "Write"
	
set output "wandering.png"
set xlabel "Threads"
set ylabel "Time (s)"
set title "Wandering details"
plot	"<awk '{if($1 == 5){print $2 \" \" $3}}' cook.dat" w lp t "Wandering function", \
 	"<awk '{if($1 == 6){print $2 \" \" $3}}' cook.dat" w lp t "Processing function"

#WINDOWS
set term x11 0
set xlabel "Threads"
set ylabel "Time (s)"
set title "Time to process"
set yrange [0:] writeback
plot 	"<awk '{if($1 == 0){print $2 \" \" $3}}' cook.dat" w lp t "Total" 
set yrange [*:*]

set term x11 1
set xlabel "Threads"
set ylabel "Time (s)"
set title "Framework pipeline"
plot 	"<awk '{if($1 == 2){print $2 \" \" $3}}' cook.dat" w lp t "Read", \
 	"<awk '{if($1 == 3){print $2 \" \" $3}}' cook.dat" w lp t "Process", \
 	"<awk '{if($1 == 4){print $2 \" \" $3}}' cook.dat" w lp t "Write"
	
set term x11 2
set xlabel "Threads"
set ylabel "Time (s)"
set title "Wandering details"
plot	"<awk '{if($1 == 5){print $2 \" \" $3}}' cook.dat" w lp t "Wandering function", \
 	"<awk '{if($1 == 6){print $2 \" \" $3}}' cook.dat" w lp t "Processing function"

pause -1
