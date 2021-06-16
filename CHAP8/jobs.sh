#!/bin/bash

# compile
echo "compiling reduced CCC code"
# gfortran -O2 -w ccc.f lapack.f iqpackd.f -o ccc


# # kgrid

# # kgrid set
# k_set="1 2"

# # theta set
# th_set="0 1 2"

# # kgrid jobs
# for k in ${k_set} ; do
#   echo "k: ${k}"
#   for th in ${th_set} ; do
#     echo "theta: ${th}"

#     dir="./data/k${k}_th${th}/"

#     echo "job: ${dir}"

#     echo "staging ccc.in"
#     cp -v "${dir}ccc.in" ./ccc.in

#     echo "executing ./ccc"
#     ./ccc < ccc.in > ccc.out

#     echo "extracting output"
#     mv -v ccc.out ${dir}
#     mv -v singlet* ${dir}
#     mv -v triplet* ${dir}
#     mv -v totalcs ${dir}

#     echo ""
#   done
# done

# book jobs

# energy set
en_set=$(echo {1..100..10})

echo ${en_set}
for en in ${en_set} ; do
  echo "energy: ${en}"

  dir="./data/book/en_${en}/"

  echo "making directory"
  mkdir -v -p ${dir}

  echo "creating ccc.in"
  sed "1 s/^[^,]*/${en}/g" book.in > "${dir}ccc.in"

  # echo "staging ccc.in"
  # cp -v "${dir}ccc.in" ./ccc.in

  # echo "executing ./ccc"
  # ./ccc < ccc.in > ccc.out

  # echo "extracting output"
  # mv -v ccc.out ${dir}
  # mv -v singlet* ${dir}
  # mv -v triplet* ${dir}
  # mv -v totalcs ${dir}

  echo ""
done
