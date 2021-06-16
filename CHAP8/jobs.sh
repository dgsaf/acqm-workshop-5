#!/bin/bash

# compile
echo "compiling reduced CCC code"
# gfortran -O2 -w ccc.f lapack.f iqpackd.f -o ccc


# # kgrid

# kgrid set
k_set="1 2"

# theta set
th_set="0 1 2"

# kgrid jobs
for k in ${k_set} ; do
  for th in ${th_set} ; do
    echo ""
    echo "k: ${k}"
    echo "theta: ${th}"

    dir="./data/kgrid/${k}_${th}/"

    # if dir already exists, job has already been run, hence skip
    if [ -d ${dir} ] ; then
      echo "already done"
      continue
    fi

    echo "making directory"
    mkdir -v -p ${dir}

    echo "creating ccc.in"
    sed "2s|[^, ]*|${th}.0|6" "kgrid${k}.in" > "${dir}ccc.in"

    echo "staging ccc.in"
    cp -v "${dir}ccc.in" ./ccc.in

    echo "executing ./ccc"
    ./ccc < ccc.in > ccc.out

    echo "extracting output"
    mv -v ccc.out ${dir}
    mv -v singlet* ${dir}
    mv -v triplet* ${dir}
    mv -v totalcs ${dir}
  done
done

# book jobs

# energy set
en_set=$(echo {1..100..1})

echo ${en_set}
for en in ${en_set} ; do
  echo ""
  echo "energy: ${en}"

  dir="./data/book/${en}/"

  # if dir already exists, job has already been run, hence skip
  if [ -d ${dir} ] ; then
    echo "already done"
    continue
  fi

  echo "making directory"
  mkdir -v -p ${dir}

  echo "creating ccc.in"
  sed "1s|[^, ]*|${en}.0|1" book.in > "${dir}ccc.in"

  echo "staging ccc.in"
  cp -v "${dir}ccc.in" ./ccc.in

  echo "executing ./ccc"
  ./ccc < ccc.in > ccc.out

  echo "extracting output"
  mv -v ccc.out ${dir}
  mv -v singlet* ${dir}
  mv -v triplet* ${dir}
  mv -v totalcs ${dir}
done
