import pandas as pd
# from Bio.Seq import Seq
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
srrs_incl = list(srr_agg[3])

libs_dict = {}

for i in range(len(srr_agg)):
    srr_acc = srr_agg[3][i]
    sub_srr_dat = srr_dat[srr_dat["SRR_acc"] == srr_acc]
    fwd_file, rev_file = sub_srr_dat["Fwd_file"][0], sub_srr_dat["Rev_file"][0]
    region, target = srr_agg[1][i], srr_agg[2][i]
    if "." in str(fwd_file) and "." in str(rev_file):
        layout = "pe"
    elif if "." in str(fwd_file):
        layout = "se"
    if target not in libs_dict and region not in libs_dict[target]:
        libs_dict[target] = {region : {srr_acc: {}}}
    else:
        libs_dict[target][region][srr_acc] = {}
    flist = glob.glob(F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/*/*{srr_acc}*')
    qiime_formatted = [x for x in flist if os.path.basename(x).count("_") == 5]
    cutadapt_formatted = [x for x in flist if os.path.basename(x).count("cutadapt") == 1]
    partially_formatted = [x for x in flist if os.path.basename(x).count("cutadapt") == 1 and os.path.basename(x).count("_") == 4]
    uncut_formatted_pe = [x for x in flist if os.path.basename(x).count("_") == 1]
    uncut_formatted_se = [x for x in flist if os.path.basename(x).count("_") == 0]
    libs_dict[target][region][srr_acc] = {"Qiime_Formatted" : qiime_formatted, "Qiime_Formatted_Count" : len(qiime_formatted),
                                          "Cutadapt_Formatted" : cutadapt_formatted, "Cutadapt_Formatted_Count" : len(cutadapt_formatted),
                                          "Partitally_Formatted" : partitally_formatted, "Partitally_Formatted_Count" : len(partitally_formatted),
                                          "Uncut_Formatted_PE" : uncut_formatted_pe, "Uncut_Formatted_PE_Count" : len(uncut_formatted_pe),
                                          "Uncut_Formatted_SE" : uncut_formatted_se, "Uncut_Formatted_SE_Count" : len(uncut_formatted_se)}
    if layout == "pe":
        libs_dict[target][region][srr_acc]["Qiime_Ready"] = (len(qiime_formatted) == 2 and qiime_formatted[0].split("/")[7] == F"{target}_{region}_{layout}")
        # libs_dict[target][region][srr_acc]["In_Correct_Location"] = (len(qiime_formatted) == 2 and cutadapt_formatted[0].split("/")[7] == F"{target}_{region}_{layout}")
    elif layout == "se":
        libs_dict[target][region][srr_acc]["Qiime_Ready"] = (len(qiime_formatted) == 1 and qiime_formatted[0].split("/")[7] == F"{target}_{region}_{layout}")


for t in libs_dict:
    for r in libs_dict[t]:
        for layout in ["pe", "se"]:
            all_seqs = glob.glob(F'/mnt/home/liberjul/He_Lab/phyllosphere_meta-analysis/SRA_files/{t}_{r}_{layout}/*.fastq.gz')
            if len(all_seqs) > 0:
                fwd = [x for x in all_seqs if "_R1_" in x]
                if "SRR" in fwd[0]:
                    acc = [x.split("SRR")[1].split("_")[0]+"SRR" for x in all_seqs if "_R1_" in x]
                elif "ERR" in fwd[0]:
                    acc = [x.split("ERR")[1].split("_")[0]+"ERR" for x in all_seqs if "_R1_" in x]
                if layout != "se":
                    rev = [x for x in all_seqs if "_R2_" in x]
                    if len(rev) == len(fwd):
                        matching_len = True
                    else:
                        matching_len = False
                else:
