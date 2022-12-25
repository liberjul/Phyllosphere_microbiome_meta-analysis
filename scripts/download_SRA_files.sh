#!/bin/bash


module load SRA-Toolkit/2.10.7-centos_linux64

FP=~/He_Lab/phyllosphere_meta-analysis/SRA_files

while IFS= read -r line
do
  TARGET=$(echo $line | cut -f1 -d' ')
  REGION=$(echo $line | cut -f2 -d' ')
  SRA_ACC=$(echo $line | cut -f3 -d' ')
  if ! [ -d $FP/"$TARGET"_"$REGION" ]
  then
    mkdir $FP/"$TARGET"_"$REGION"
  fi
  echo "$TARGET"_"$REGION/$SRA_ACC"
  cd $FP/"$TARGET"_"$REGION"
  fasterq-dump $SRA_ACC
done < $1

cd $FP/../scripts
