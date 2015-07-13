reset

mean(x,y) = ((x*y+100*(16384-y))/16384)
serial(x) = 16384*100/(100*x+1*(16384-x))
derak(x)  = x == 6400?x:-1000000

WEIGHT = 2

set terminal png enhanced font "Arial,12" dashed size 800, 600

set output "Nxy_i_10.png"
set multiplot
set origin 0.01,0.01
set size 0.99, 0.99
set xrange [ 0 : 1 ]
set yrange [ 0.01 : 0.33 ]
set xtics 0,0.1,1#0,1024,16384
set ytics 0.01,0.01,1
set key spacing 0.0#2.4
set key top left
# \n Рис.3. Сравнение макротеплопроводности ячейки для \n включений различной формы
set xlabel "Коэффициент армирования {с}" font "Monospace, 30"
set ylabel "Коэффициент Пуассона {/Symbol n_{/Nosymbol xy}}" font "Arial, 30"# offset 2

set style line 1 lt 1 pt 1 lw 2 lc rgb "black"
set style line 2 lt 2 pt 1 lw 2 lc rgb "black" 
set style line 3 lt 3 pt 1 lw 2 lc rgb "black" 
set style line 4 lt 4 pt 1 lw 2 lc rgb "black" 
set style line 5 lt 1 pt 1 lw 2 lc rgb "red" 
set style line 6 lt 1 pt 2 lw 2 lc rgb "red" 
set style line 7 lt 1 pt 1 lw 2 lc rgb "red" 
set style line 8 lt 1 pt 2 lw 2 lc rgb "red"  

plot \
 "../concrete/mata-quadrate.gpd" using ($1/16384):($3) w l ls 1 ti "{/Arial=10 квадрат}",\
 "../concrete/mata-shell.gpd" using ($1/16384):($3) w l ls 2 ti "{/Arial=10 трубка}",\
 "../concrete/mata-cross.gpd" using ($1/16384):($3) w l ls 3 ti "{/Arial=10 крестовина}", \
 "../concrete/mata-circ.gpd" using ($1/16384):($3) w l ls 4 ti "{/Arial=10 круг}", \
 "../../ch_fork.gpd" using ($1/16384):($3) w lp ls 5 ti "{/Arial=10 нижня оценка}" , \
 "../../ch_fork.gpd" using ($1/16384):($15) w lp ls 6 ti "{/Arial=10 верхняя оценка}", \
 "mata-hexagon.gpd" using ($1):($3) w lp lt 2 ti "{/Arial=10 hex}", \
 "mata-hexagon.gpd" using ($1):($14) w lp lt 3 ti "{/Arial=10 hex}", \
 "../../ch_fork.gpd" using ($1/16384):(($13*$10 - $11*$10 - 4.0 * $8**2 * $11 * $13) / ($13*$10 + $11*$10 + 4.0 * $8**2 * $11 * $13)) w lp ls 6 ti "{/Arial=10 верхняя оценка}"

#, "all/l" with labels font "Monospace, 20" ti ""#, "plas/unconst area/7" using (derak($1)/16384):($2/100) with points ti "" lt -1 pt 1 ps 3#, "shell/unconst area/7" using (derak($1)/16384):($2/100) with points ti "{/Monospace=20 S_в = 0.39}" lt -1 pt 1 ps 3, "all/l1" with labels font "Monospace, 20" ti ""

#set titl "Относительное отклонение от среднего \n ({/Symbol=20 l_s} - {{/Symbol=20 l}_с}) / {{/Symbol=20 l}_с} " font "Monospace, 20"
#set origin 0.54,0.43
#set size 0.4, 0.4
#set noxlabel 
#set xlabel "{S_B/S}" 
#set noylabel
#set xrange [ 0 : 1 ]
#set yrange [ 0 : 1 ]
#set xtics 0,0.1,1
#set ytics 0,0.1,1
#set grid x y
#plot "cube/unconst area/6" using ($1/16384):(((mean(1,$1)-$2))/mean(1,$1)) with lines ti "" lw 3, "shell/unconst area/7" using ($1/16384):(((mean(1,$1)-$2))/mean(1,$1)) with lines ti "" lw 3, "plas/unconst area/7" using ($1/16384):(((mean(1,$1)-$2))/mean(1,$1)) with lines ti "" lw 3







