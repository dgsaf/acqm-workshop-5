#
set terminal epslatex input color solid


# files
dir(n)=sprintf("../data/output/1s_%is.dir.txt", n)
exc(n)=sprintf("../data/output/1s_%is.exc.txt", n)

on_dir(n)=sprintf("../data/output/1s_%is.dir.on.txt", n)
on_exc(n)=sprintf("../data/output/1s_%is.exc.on.txt", n)


# common settings
set palette defined (0 "blue" , 1 "red")
unset colorbox
set grid xtics ytics ztics
set key \
  top right \
  box opaque \
  samplen 1 spacing 0.6 height +0.6


# surface plot settings
set xlabel "$k_{i}$"
set ylabel "$k_{f}$"
set zlabel "$V(k_{f}, k_{i})$"

set xrange [0:5]
set yrange [0:5]

set format x "\\scriptsize %.0f"
set format y "\\scriptsize %.0f"
set format z "\\scriptsize %.2f"

set view 60,120

unset key


# figure: direct matrix element surface plot for each (1s -> ns) transition
do for [n=1:3] {
   set output sprintf("figure_1s_%is_dir.tex", n)


   set title sprintf("Direct Matrix Elements [$1s \\to %is$]", n)

   set autoscale z

   splot dir(n) using 1:2:3 title "Direct" with lines lc "red"

   set output
}

# figure: exchange matrix element surface plot for each (1s -> ns) transition
do for [n=1:3] {
   set output sprintf("figure_1s_%is_exc.tex", n)

   # set key width -3.75

   set title sprintf("Exchange Matrix Elements [$1s \\to %is$]", n)

   splot exc(n) using 1:2:3 title "Exchange" with lines lc "blue"

   set output
}