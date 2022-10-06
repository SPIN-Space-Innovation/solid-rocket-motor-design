f(x) = 1
plot f(x), 'exit_pressure_NSE.dat' with line, 'exit_pressure_sup.dat' with line
set grid
set xlabel 'time [sec]'
set ylabel 'pressure [bar]'
set title 'pressure'
set term png size 1500,1000
set output 'exitPressure.png'
replot

plot 'ChamberPressure.dat' with line
set grid
set xlabel 'time [sec]'
set ylabel 'pressure [Pa]'
set title 'pressure'
set term png size 1500,1000
set output 'ChamberPressure.png'
replot


plot 'massFlowRate.dat' with line
set grid
set xlabel 'time [sec]'
set ylabel 'm_dot [kg/sec]'
set title 'm_dot'
set term png size 1500,1000
set output 'm_dot.png'
replot

plot 'exitTemperature.dat' with line
set grid
set xlabel 'time [sec]'
set ylabel 'T_e [K]'
set title 'm_dot'
set term png size 1500,1000
set output 'T_e.png'
replot

plot 'thrust.dat' with line
set grid
set xlabel 'time [sec]'
set ylabel 'thrust [N]'
set title 'thrust'
set term png size 1500,1000
set output 'thrust.png'
replot