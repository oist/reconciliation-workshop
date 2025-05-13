#!/bin/bash
commit1=$!
commit2=$2
cores=40

if [ "$#" -ne 2 ]; then
  echo "Syntax: "
  echo "compare.sh commt1 commit2"
  exit 2
fi

if [ "$1" = "$2" ]; then
  echo " First and second argument must be different"
  exit 3
fi

giturl="https://codeberg.org/Exelixis-Lab/coraxlib.git"

echo "Comparing commits:"
echo "- $1"
echo "- $2"



temp="temp_bench"
rm -rf $temp
mkdir $temp
cd $temp

commit1=$1
commit2=$2
repo1="$commit1"
repo2="$commit2"

buildrepo() {
  commit=$1
  repo=$commit
  out="out_$commit.json"
  echo $commit
  git clone --recursive  $giturl $repo || { echo "Cannot clone branch $commit" ; exit 1; }
  cd $repo
  git checkout $commit
  mkdir build
  cd build
  cmake -DCORAX_TESTS_ENABLED=OFF -DCORAX_BENCH_ENABLED=ON ..
  make -j $cores
  cd ..
  ./build/bench/src/corax-bench --benchmark_out=$out
  cd ..

}


buildrepo $commit1
buildrepo $commit2

out1=$commit1/out_$commit1.json
out2=$commit2/out_$commit2.json

../bench/gbench_tools/compare.py benchmarks $out1 $out2
echo  $out1

