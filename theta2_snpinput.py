#__author__ = "Janice"

import os
import re
import pandas as pd

newdir = "C:\Users\Owner\Box Sync\\theta2\\"
snpdir = "C:\Users\Owner\Box Sync\mutect_output\mutect_output_raw\\"


for file in os.listdir(snpdir):
    if file.endswith("callstats.txt"):

        snpmutect=pd.read_table(snpdir + file, header=1, sep="\t")
        #snpfilter=snpmutect[(snpmutect.judgement=="KEEP")]
        snpfilter=snpmutect
        #remove 'chr' from chromosome names
        snpfilter['contig']=snpfilter['contig'].str[3:]

        #remove random chrN_random tables they mess up theta
        boolist=snpfilter.contig.str.contains('(_random)')
        notboolist=[not i for i in boolist]
        snpclean=snpfilter[notboolist]

        #add a fucking hashtag to the first line, because god and theta2 hate me
        snpclean=snpclean.rename(columns = {'contig':'#contig'})
        tumor=snpclean[['#contig','position','t_ref_count','t_alt_count']]
        normal=snpclean[['#contig','position','n_ref_count','n_alt_count']]

        tumor.to_csv(newdir +  re.sub(r"dedup_realigned_recalib_callstats.txt$|merged_sort_callstats.txt$", 'tumor_SNP.txt', file), sep="\t", index=False)
        normal.to_csv(newdir + re.sub(r"dedup_realigned_recalib_callstats.txt$|merged_sort_callstats.txt$", 'norm_SNP.txt', file), sep="\t", index=False)
