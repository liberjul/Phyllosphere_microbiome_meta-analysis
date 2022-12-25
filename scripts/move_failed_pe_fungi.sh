#!/bin/bash --login

FP=~/He_Lab/phyllosphere_meta-analysis/SRA_files

while IFS= read -r line
do
  if [ $line == *"Fungi"* ]; then
    TARGET=$(echo $line | cut -f1 -d' ')
    REGION=$(echo $line | cut -f2 -d' ')
    SAMPLE=$(echo $line | cut -f3 -d' ')
    mkdir -p $FP/"$TARGET"_"$REGION"_se
    cp $FP/"$TARGET"_"$REGION"_pe/$SAMPLE*R1* $FP/"$TARGET"_"$REGION"_se/
  fi
done < no_reads_fungi.txt
