#/usr/bash

gnuplot --persist -e 'plot "data.txt" using 2:1 title "rotation" with lines,"data.txt" using 3:1 title "vertical disp" with lines'
#gnuplot --persist -e 'plot "data.txt" using 4'
#gnuplot --persist -e 'plot "data.txt" using 5'
