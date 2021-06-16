set terminal epslatex input color solid

# files
triplet(k,th)=sprintf("../data/kgrid/%i_%i/triplet.11", k, th)

# common settings
set palette defined (0 "blue" , 1 "red")
unset colorbox
set grid xtics ytics ztics
set key \
  top right \
  box opaque \
  samplen 1 spacing 0.6 height +0.6

# 2d plot settings
set xlabel "$k$"
set ylabel "Matrix Elements"

set xrange [0:*]
set yrange [-0.5:0.5]

set format x "\\scriptsize %.0f"
set format y "\\scriptsize %.2f"
set key width -3.5 spacing 0.8 height +0.6

# figure(s): on-shell matrix element plot for each (1s -> ns) transition
do for [th = 0:2] {
  set output sprintf("ccc_figure_triplet_%i.tex", th)

  set title \
    sprintf("%s [$S = 1$, $1s \\to 1s$, $\\theta = %i$]"\
    , "Half-on-Shell Matrix Elements", th)

  plot \
    for [k=1:2] \
      triplet(k,th) using 1:2 smooth unique \
      title sprintf("\\scriptsize $K(k', k)$, k-grid %i", k) \
      with lines palette frac (0.0+((k-1.0)/1.0)) , \
    for [k=1:2] \
      triplet(k,th) using 1:3 smooth unique \
      title sprintf("\\scriptsize $V(k', k)$, k-grid %i", k) \
      with lines palette frac (0.0+((k-1.0)/1.0)) dashtype 2

  set output
}
