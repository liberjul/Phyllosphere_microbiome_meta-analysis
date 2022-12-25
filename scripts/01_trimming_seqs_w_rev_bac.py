import pandas as pd
from Bio.Seq import Seq
import subprocess, os, glob

primer_dict = {
"338F-806R" : ["ACTCCTACGGGAGGCAGCA", "GGACTACNVGGGTWTCTAAT"],
"341F-785R" : ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"],
"515F-806R" : ["GTGYCAGCMGCCGCGGTAA", "GGACTACNVGGGTWTCTAAT"],
"341F-805R" : ["CCTACGGGNGGCWGCAG", "GACTACHVGGGTATCTAATCC"],
"5.8SFun-ITS4Fun" : ["AACTTTYRRCAAYGGATCWCT", "AGCCTCCGCTTATTGATATGCTTAART"],
"501F-706R" : ["CAGCCTACGGGNGGCWGCAG", "ACTACHVGGGTATCTAATCC"],
"515F-926R" : ["GTGYCAGCMGCCGCGGTAA", "CCGYCAATTYMTTTRAGTTT"],
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
srr_agg = pd.read_csv("srr_aggregatted_by_target_region.txt", sep = "\t", header = None)

for srr_acc in list(srr_agg[3]):
    sub_srr_dat = srr_dat[srr_dat["SRR_acc"] == srr_acc].to_dict()
    index = list(sub_srr_dat["SRR_acc"].keys())[0]
    fwd_file, rev_file = sub_srr_dat["Fwd_file"][index], sub_srr_dat["Rev_file"][index]
    target, region = sub_srr_dat["Target"][index], sub_srr_dat["Region"][index]
    primer_set = sub_srr_dat["Primer_set"][index]
    fwd_primer_seq, rev_primer_seq = primer_dict[primer_set]
    fwd_rc, rev_rc = str(Seq(fwd_primer_seq).reverse_complement()), str(Seq(rev_primer_seq).reverse_complement())
    if target == "Bacteria":
        if "." in str(fwd_file) and "." in str(rev_file):
            inpath = F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/uncut_seqs/{srr_acc}'
            outpath = F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/{target}_{region}_pe/JL_{srr_acc}_L001_R'
            if len(glob.glob(F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/{target}_{region}_pe/JL_...._{srr_acc}_L001_R._001.fastq.gz')) > 0:
                pass
            elif os.path.exists(inpath + "_1.sampled.fastq.gz"):
                if primer_set in ["926F-1392R", "799F-1492R"] or region == "v3-v4" or region == "ITS1-ITS2":
                    cmd = F"cutadapt -g {fwd_primer_seq} -G {rev_primer_seq} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz -p {outpath}2_001.cutadapt.fastq.gz {inpath}_1.sampled.fastq.gz {inpath}_2.sampled.fastq.gz"
                    subprocess.run(cmd, shell = True)
                else:
                    cmd = F"cutadapt -a {fwd_primer_seq}...{rev_rc} -A {rev_primer_seq}...{fwd_rc} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz -p {outpath}2_001.cutadapt.fastq.gz {inpath}_1.sampled.fastq.gz {inpath}_2.sampled.fastq.gz"
                    subprocess.run(cmd, shell = True)
            elif os.path.exists(inpath + "_1.fastq.gz") and os.path.exists(inpath + "_2.fastq.gz") and not os.path.exists(F"{outpath}1_001.cutadapt.fastq.gz"):
                if primer_set in ["926F-1392R", "799F-1492R"] or region == "v3-v4" or region == "ITS1-ITS2":
                    cmd = F"cutadapt -g {fwd_primer_seq} -G {rev_primer_seq} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz -p {outpath}2_001.cutadapt.fastq.gz {inpath}_1.fastq.gz {inpath}_2.fastq.gz"
                    subprocess.run(cmd, shell = True)
                else:
                    cmd = F"cutadapt -a {fwd_primer_seq}...{rev_rc} -A {rev_primer_seq}...{fwd_rc} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz -p {outpath}2_001.cutadapt.fastq.gz {inpath}_1.fastq.gz {inpath}_2.fastq.gz"
                    subprocess.run(cmd, shell = True)
        elif "." in str(fwd_file):
            inpath = F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/{srr_acc}'
            outpath = F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/{target}_{region}_se/JL_{srr_acc}_L001_R'
            if os.path.exists(inpath + ".fastq.gz") and not os.path.exists(F"{outpath}1_001.cutadapt.fastq.gz"):
                cmd = F"cutadapt -a {fwd_primer_seq}...{rev_rc} -e .1 -j 24 -m 150 -o {outpath}1_001.cutadapt.fastq.gz {inpath}.fastq.gz"
                subprocess.run(cmd, shell = True)
        else:
            print(F"{SRR_acc} does not have a filename")
