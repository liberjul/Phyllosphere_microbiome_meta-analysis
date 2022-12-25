#!/bin/bash --login

FP=~/He_Lab/phyllosphere_meta-analysis/SRA_files

TARGET="Bacteria"
for REGION in "v1-v3" "v3-v4" "v4" "v4-v5" "v5-v6" "v6-v7"
do
  for LAYOUT in "se" "pe"
  do
    i="$TARGET"_"$REGION"_"$LAYOUT"
    if [ -d $FP/"$i" ]; then
      sed "s/<target>/$TARGET/" 04_classify_constax.sb | sed "s/<region>/$REGION/" | sed "s/<layout>/$LAYOUT/" > 04_classify_constax_"$i".sb
      echo "$i"
      # sbatch 03_qiime_import_dada2_denoise_"$i".sb
    fi
  done
done
TARGET="Fungi"
for REGION in "ITS1" "ITS2" "ITS1-ITS2"
do
  for LAYOUT in "se" "pe"
  do
    i="$TARGET"_"$REGION"_"$LAYOUT"
    if [ -d $FP/"$i" ]; then
      sed "s/<target>/$TARGET/" 04_classify_constax.sb | sed "s/<region>/$REGION/" | sed "s/<layout>/$LAYOUT/" > 04_classify_constax_"$i".sb
      echo "$i"
      # sbatch 03_qiime_import_dada2_denoise_"$i".sb
    fi
  done
done
