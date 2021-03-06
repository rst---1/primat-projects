reset

set terminal png enhanced font "Arial,12" dashed dl 5 size 800, 600

set output "Теплопроводность вдоль оси Ox до 035.png"
set origin 0.01, 0.01
set size 0.99, 0.99
set xrange [ 0.00 : 0.35 ]
set key spacing 2.4
#set nokey
set key top left 
set key samplen 14

set style line 1 lt 1 pt 1 lw 2 lc rgb "red"
set style line 2 lt 2 pt 1 lw 2 lc rgb "green" 
set style line 3 lt 3 pt 1 lw 2 lc rgb "blue" 
set style line 4 lt 1 pt 4 ps 2.0 lw 2 lc rgb "black" 
set style line 5 lt 1 pt 1 ps 2.0 lw 2 lc rgb "black" 
set style line 6 lt 1 pt 2 ps 2.0 lw 2 lc rgb "black" 

set xlabel "Коэффициент армирования" font "Monospace, 20"
set ylabel "Коэффициент теплопроводности {{/Symbol=20 l}_{xx}}" font "Monospace, 20" offset 2

plot \
 "meta-hex_before_035.gpd" u 1:2 w l ti "1" lw 2 lt 1, \
 "jancovsky_before_035.gpd" u 1:2 w l ti "2" lw 2 lt 2, \
 "jancovsky_before_035.gpd" u 1:4 w l ti "3" lw 2 lt 3
 

set output "Теплопроводность вдоль оси Oy до 035.png"
set ylabel "Коэффициент теплопроводности {{/Symbol=20 l}_{yy}}" font "Monospace, 20" offset 2

plot \
 "meta-hex_before_035.gpd" u 1:3 w l ti "1" lw 2 lt 1, \
 "jancovsky_before_035.gpd" u 1:3 w l ti "2" lw 2 lt 2, \
 "jancovsky_before_035.gpd" u 1:5 w l ti "3" lw 2 lt 3
