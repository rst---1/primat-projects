reset

WEIGHT = 2
DEVIDER = 100

arithmetic(s,l1,l2) = ((s*l2+l1*(1.0-s))/1.0)
harmonic(s,l1,l2) = 1.0 / (s/l2 + (1.0 - s)/l1)
serial(x) = 16384*DEVIDER/(DEVIDER*x+1*(16384-x))
derak(x)  = x == 6400?x:-1000000

set terminal png enhanced font "Serif,12" size 800, 600

set output "4.png"
set multiplot
set origin 0.00,0.00
set size 1.0, 1.0
set xrange [ 0 : 1 ]
set yrange [ 0.01 : 1 ]
set xtics 0,0.1,1#0,1024,16384
set ytics 0.01,0.11,1
set key spacing 2.1
set key top right
# set title "Sine [-Pi..Pi]"
# Отношение площади включения к общей площади {S_B/S}
set xlabel "Коэффициент армирования {/Symbol \161^{/Nosymbol I}}" font "Serif, 30"
set ylabel "Коэффициент теплопроводности {{/Symbol l}_{xx}}" font "Serif, 20"

set style line 1 lt 1 pt 1 lw 2 lc rgb "#000000"
set style line 2 lt 2 pt 1 lw 2 lc rgb "#00FFFF" 
set style line 3 lt 3 pt 1 lw 2 lc rgb "#FF00FF" 
set style line 4 lt 4 pt 1 lw 2 lc rgb "#FF7F00" 
set style line 5 lt 1 pt 5 ps 1.5 lw 2 lc rgb "red" 
set style line 6 lt 1 pt 7 ps 1.5 lw 2 lc rgb "red" 
set style line 7 lt 1 pt 10 lw 2 lc rgb "blue" 


#set titl "Относительное отклонение \n от круглого включения" font "Serif, 15"
set origin 0.43,0.40
set size 0.37, 0.42
set noxlabel 
set xlabel "Гексагональная упаковка"  font "Serif, 12"
set noylabel
set xrange [ 0 : 1 ]
set yrange [ 0 : 1.0 ]
set xtics 0,0.2,1
set ytics 0,0.2,1.0
set grid x y
plot \
 "../mata-hexagon.gpd" u ($1):(($5 - $2) / $2) w l ls 2 ti "" , \
 "../mata-hexagon.gpd" u ($1):(($6 - $2) / $2) w l ls 3 ti ""

#set titl "Относительное отклонение \n от круглого включения" font "Serif, 15"
set origin 0.08,0.40
set size 0.37, 0.42
set noxlabel 
set xlabel "Тетрагональная упаковка"  font "Serif, 12"
set noylabel
set xrange [ 0 : 1 ]
set yrange [ 0 : 1.0 ]
set xtics 0,0.2,1
set ytics 0,0.2,1.0
set grid x y
plot \
 "../mata-circ.gpd" u ($1):(($5 - $2) / $2) w l ls 2 ti "" , \
 "../mata-circ.gpd" u ($1):(($6 - $2) / $2) w l ls 3 ti ""





