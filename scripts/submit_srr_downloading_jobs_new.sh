#!/bin/bash --login

while IFS= read -r line
do
  TARGET=$(echo $line | cut -f1 -d' ')
  REGION=$(echo $line | cut -f2 -d' ')
  sed "s/<target>/$TARGET/" download_sra_files_job_new.sb | sed "s/<region>/$REGION/"> download_sra_files_job_"$TARGET"_"$REGION"_new.sb
  grep "$TARGET.$REGION" srr_aggregatted_by_target_region_new_downloads.txt > srr_list_"$TARGET"_"$REGION"_new.txt
  # sbatch download_sra_files_job_"$TARGET"_"$REGION"_new.sb
done < srr_targets_regions.txt
