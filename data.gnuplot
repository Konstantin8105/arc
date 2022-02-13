#/usr/bash

gnuplot --persist -e 'plot "data.txt" using 2:1 title "rotation" with lines,"data.txt" using 3:1 title "vertical disp" with lines'
#gnuplot --persist -e 'plot "data.txt" using 4'
#gnuplot --persist -e 'plot "data.txt" using 5'
#
gnuplot --persist -e 'plot "arc3.txt" using 2:1 title "1d" with lines, "arc3.txt" using 2:3 title "integral" with line, "arc3.txt" using 2:4 title "perfect" with line'
gnuplot --persist -e 'plot "arc4.txt" using 2:1 title "u" with lines, "arc4.txt" using 3:1 title "v" with lines'
