#/usr/bash

gnuplot --persist -e 'plot "data.txt" using 2:1 ,"data.txt" using 3:1'
gnuplot --persist -e 'plot "data.txt" using 4'
gnuplot --persist -e 'plot "data.txt" using 5'
