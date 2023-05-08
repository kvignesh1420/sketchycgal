#!/bin/bash

# function to download and extract tar.gz file
download_extract() {
  url=$1
  filename=$(basename "$url")
  dirname="${filename%.*.*}"
  mkdir -p "$2/$dirname"
  wget -q --show-progress "$url" -P "$2/$dirname"
  tar -xzf "$2/$dirname/$filename" -C "$2/$dirname"
  rm "$2/$dirname/$filename"
}

# download G67 data files for experiments

for i in {1..67}; do
  url="https://www.cise.ufl.edu/research/sparse/MM/Gset/G$i.tar.gz"
  download_extract "$url" "data/G"
done

# selectively download large scale graphs from SuiteSparse collection
# NOTE: SuiteSparse is included in the DIMACS10 collection

url=https://suitesparse-collection-website.herokuapp.com/MM/DIMACS10/delaunay_n20.tar.gz
download_extract "$url" "data/DIMACS10"