#
set terminal epslatex input color solid


# set energy (manually change this for each energy in {0.1, 10.0})
en=10.0


# files
dir(n)=sprintf("../data/output/%.1f_1s_%is.dir.txt", en, n)
exc(n)=sprintf("../data/output/%.1f_1s_%is.exc.txt", en, n)

dir_on(n)=sprintf("../data/output/%.1f_1s_%is.dir.on.txt", en, n)
exc_on(n)=sprintf("../data/output/%.1f_1s_%is.exc.on.txt", en, n)


# common settings
set palette defined (0 "blue" , 1 "red")
unset colorbox
set grid xtics ytics ztics
set key \
  top right \
  box opaque \
  samplen 1 spacing 0.6 height +0.6


# surface plot settings
set xlabel "$k'$"
set ylabel "$k$"
set zlabel "$V(k', k)$"

set xrange [0:5]
set yrange [0:5]
set autoscale z

set format x "\\scriptsize %.0f"
set format y "\\scriptsize %.0f"
set format z "\\scriptsize %.2f"

set view 60,120
unset key


# figure(s): direct matrix element surface plot for each (1s -> ns) transition
do for [n=1:3] {
  set output sprintf("figure_%.0f_1s_%is_dir.tex", en, n)

  set title sprintf("$D_{%i, 1}(k', k)$, $E = \\SI{%.0f}{\\eV}$", n, en)

  splot dir(n) using 1:2:3 title "Direct" with lines lc "red"

  set output
}


# figure(s): exchange matrix element surface plot for each (1s -> ns) transition
do for [n=1:3] {
  set output sprintf("figure_%.0f_1s_%is_exc.tex", en, n)

  # set title sprintf("Exchange Matrix Elements [$1s \\to %is, E = %.0f$]", n, en)
  set title sprintf("$X_{%i, 1}(k', k)$, $E = \\SI{%.0f}{\\eV}$", n, en)

  splot exc(n) using 1:2:3 title "Exchange" with lines lc "blue"

  set output
}


# 2d plot settings
set xlabel "$k$"
set ylabel "$V(k', k)$"

set xrange [0:5]
set autoscale y

set format x "\\scriptsize %.0f"
set format y "\\scriptsize %.2f"

set key width -3.5 spacing 0.8 height +0.6

# figure(s): on-shell matrix element plot for each (1s -> ns) transition
do for [n=1:3] {
  set output sprintf("figure_%.0f_1s_%is_on.tex", en, n)

  set title sprintf("On-Shell Matrix Elements [$1s \\to %is$]", n)

  plot \
    dir_on(n) using 1:2 \
      title sprintf("\\scriptsize $D_{%i, 1}(k', k)$", n) \
      with lines lc "red", \
    dir_on(n) using 1:3 \
      title sprintf("\\scriptsize [Analytic] $D_{%i, 1}(k', k)$", n) \
      with lines lc "black", \
    exc_on(n) using 1:2 \
      title sprintf("\\scriptsize $X_{%i, 1}(k', k)$", n) \
      with lines lc "blue", \
    exc_on(n) using 1:3 \
      title sprintf("\\scriptsize [Analytic] $X_{%i, 1}(k', k)$", n) \
      with lines lc "black"

  set output
}
