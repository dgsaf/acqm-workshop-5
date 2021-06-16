#!/bin/bash

# kgrid set
k_set="1 2"

# theta set
th_set="0 1 2"

# compile
gfortran -O2 -w ccc.f lapack.f iqpackd.f -o ccc

for k in ${kset}; do
  for th in ${th_set}; do
    dir="./data/k${k}_th${th}/"

    cp "${dir}ccc.in" ./ccc.in

    ./ccc < ccc.in > ccc.out

    mv ccc.out "${dir}"
    mv singlet* "${dir}"
    mv triplet* "${dir}"
    mv totalcs "${dir}"
  done
done
