f(x) = 1
plot f(x) title 'ambient pressure', 'exit_pressure_sup.dat' with line, 'exit_pressure_NSE.dat' with line
set grid
set xlabel 'time [sec]'
set ylabel 'pressure [bar]'
set title 'pressure'
set term png size 1500,1000
set output 'exit_pressure_2.png'
replot
