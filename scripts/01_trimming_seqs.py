import pandas as pd
import subprocess, os

primer_dict = {
"338F-806R" : ["ACTCCTACGGGAGGCAGCA", "GGACTACNVGGGTWTCTAAT"],
"341F-785R" : ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"],
"515F-806R" : ["GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"],
"341F-805R" : ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"],
"5.8SFun-ITS4Fun" : ["AACTTTYRRCAAYGGATCWCT", "AGCCTCCGCTTATTGATATGCTTAART"],
"501F-706R" : ["TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCCTACGGGNGGCWGCAG", "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGACTACHVGGGTATCTAATCC"],
"515F-926R" : ["GTGYCAGCMGCCGCGGTAA", "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGATCCGYCAATTYMTTTRAGTTT"],
"785F-1064R" : ["GGATTAGATACCC", "CGACRRCCATGCANCACCT"],
"799F-1115R" : ["AACMGGATTAGATACCCKG", "AGGGTTGCGCTCGTTG"],
"799F-1193R" : ["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"],
"799F-1492R" : ["AACMGGATTAGATACCCKG", "TACCTTGTTACGACTT"],
"926F-1392R" : ["AAACTYAAAKGAATTGACGG", "ACGGGCGGTGTGTRC"],
"B341F-B806R" : ["CCTACGGGAGGCAGCAG", "GGACTACHVGGGTWTCTAAT"],
"B799F-B1194R" : ["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"],
"BAC799F-BAC1193R" : ["AACMGGATTAGATACCCKG", "ACGTCATCCCCACCTTCC"],
"fITS7-ITS4" : ["GTGAATCATCGAATCTTTG" , "TCCTCCGCTTATTGATATGC"],
"fITS7-ITS4R" : ["GTGAATCATCGAATCTTTG", "TCCTCCGCTTATTGATATGC"],
"Illumina V1-V3" : ["AGAGTTTGATCCTGGCTCAG", "ATTACCGCGGCTGCTGG"],
"ITS-1F-F-ITS1-1F-R" : ["CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"],
"ITS1F-ITS2" : ["CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"],
"ITS1F-ITS2R" : ["CTTGGTCATTTAGAGGAAGTAA", "GCTGCGTTCTTCATCGATGC"],
"ITS1f-ITS4" : ["TTGGTCATTTAGAGGAAGTA", "TCCTCCGCTTATTGATATGC"],
"ITS1F-KYO1-ITS2-KYO1" : ["CTHGGTCATTTAGAGGAASTAA", "TAGAGGAAGTAAAAGTCGTAA"],
"ITS1O-5.8s-O-R" : ["CGGAAGGATCATTACCAC", "AGCCTAGACATCCACTGCTG"],
"ITS2-F-ITS2-R" : ["GCATCGATGAAGAACGC", "CCTCCGCTTATTGATATGC"],
"ITS7-ITS4" : ["GTGAATCATCGAATCTTTG", "TCCTCCGCTTATTGATATGC"],
"TS1F12-5.8S" : ["GAACCWGCGGARGGATCA", "CGCTGCGTTCTTCATCG"],
}

srr_dat = pd.read_csv("srr_with_host_taxonomy.txt", sep = "\t")

for i in range(len(srr_dat)):
    primer_set = srr_dat["Primer_set"][i]
    fwd_primer_seq, rev_primer_seq = primer_dict[primer_set]
    fwd_file, rev_file = srr_dat["Fwd_file"][i], srr_dat["Rev_file"][i]
    # print(fwd_file, rev_file)
    if "." in str(fwd_file) and "." in str(rev_file):
        path = F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/{srr_dat["Target"][i]}_{srr_dat["Region"][i]}/{srr_dat["SRR_acc"][i]}'
        if os.path.exists(path + "_1.fastq"):
            cmd1 = F"cutadapt -g {fwd_primer_seq} -e .1 -o {path}_1.cutadapt.fastq {path}_1.fastq -j 24"
            subprocess.run(cmd1, shell = True)
        if os.path.exists(path + "_2.fastq"):
            cmd2 = F"cutadapt -g {rev_primer_seq} -e .1 -o {path}_2.cutadapt.fastq {path}_2.fastq -j 24"
            subprocess.run(cmd2, shell = True)
    elif "." in str(fwd_file):
        path = F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/{srr_dat["Target"][i]}_{srr_dat["Region"][i]}/{srr_dat["SRR_acc"][i]}'
        if os.path.exists(path + ".fastq"):
            cmd = F"cutadapt -g {fwd_primer_seq}.. -e .1 -o {path}.cutadapt.fastq {path}.fastq -j 24"
            subprocess.run(cmd, shell = True)
    else:
        print(F"{SRR_acc} does not have a filename")
