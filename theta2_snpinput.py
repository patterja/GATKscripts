#__author__ = "Janice"

import os
import re
import pandas as pd

newdir = "C:\Users\Owner\Box Sync\\theta2\\"
snpdir = "C:\Users\Owner\Box Sync\mutect_output\mutect_output_raw\\"


for file in os.listdir(snpdir):
    if file.endswith("callstats.txt"):
        snpmutect=pd.read_table(snpdir + file, header=1, sep="\t")
        snpfilter=snpmutect[(snpmutect.judgement=="KEEP")]
        tumor=snpfilter[['contig','position','t_ref_count','t_alt_count']]
        normal=snpfilter[['contig','position','n_ref_count','n_alt_count']]

        tumor.to_csv(newdir +  re.sub(r"dedup_realigned_recalib_callstats.txt$|merged_sort_callstats.txt$", 'tumor_SNP.txt', file), sep="\t", index=False)
        normal.to_csv(newdir + re.sub(r"dedup_realigned_recalib_callstats.txt$|merged_sort_callstats.txt$", 'norm_SNP.txt', file), sep="\t", index=False)
