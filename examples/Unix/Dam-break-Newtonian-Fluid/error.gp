set datafile separator ','
set logscale y 10
set key autotitle columnhead
set xlabel 'Steps'
set ylabel 'Total Error / Relative Error'

set y2tics # enable second axis
set ytics nomirror # dont show the tics on that side
set y2label "Number of iterations" # label for second axis
set y2range [0:10]

while (1) {
    plot './Resultados/Stats_Solver.csv' using 2 with lines , \
    './Resultados/Stats_Solver.csv' using 3 with lines , \
    './Resultados/Stats_Solver.csv' using 1 with lines axis x1y2 lt -1
    pause 1      # waiting time in seconds
}
