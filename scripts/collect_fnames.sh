#!/bin/bash

FP=~/He_Lab/phyllosphere_meta-analysis/SRA_files

for FOLDER in "Bacteria_v1-v3" "Bacteria_v3-v4" "Bacteria_v4" "Bacteria_v4-v5" "Bacteria_v5-v6" "Bacteria_v6-v7" "Fungi_ITS1-ITS2" "Fungi_ITS2" "Fungi_ITS1"
do
  for LAYOUT in "pe" "se"; do
    ll $FP/"$FOLDER"_"$LAYOUT"/ > $FP/../scripts/fnames_"$FOLDER"_"$LAYOUT".txt
  done
done

for FOLDER in "Bacteria_v4" "Bacteria_v5-v6"
do
  for LAYOUT in "pe" "se"; do
    for NUM in {1..5}; do
      ll $FP/"$FOLDER"_"$LAYOUT"_"$NUM"/ > $FP/../scripts/fnames_"$FOLDER"_"$LAYOUT"_"$NUM".txt
    done
  done
done
