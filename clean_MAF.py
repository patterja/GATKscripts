#!/usr/bin/env python
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Output MAF file from Oncotator output.
Usage: python clean_MAF.py
#
#Modules can be used, not used because Indel and SNP mutect output are completely different
#because I'm cursed. Something fixable, when I'm smarter.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import os
import re
import pandas as pd

newdir = "C:\Users\Owner\Box Sync\mutsig\inputs\\"
snpdir = "C:\Users\Owner\Box Sync\oncotator\oncotator_snps\\"
indeldir = "C:\Users\Owner\Box Sync\oncotator\oncotator_indel\\"

for subdir, dirs, files in os.walk(snpdir):
    #input the directory you want to walk through full of the oncotate files this is the raw everything file
    fh = open(snpdir + files[1], 'r')
    lines = fh.readlines()
    # write column names/header into the master file
    fw = open(newdir + "gistsnp_somatic_oncotate_unlabel.txt", 'w')
    fw.write(str(lines[100]))
    fh.close()
    norm_match=[]
    tumor_samp=[]
    for file in files:
        fh = open(snpdir + file, 'r')  # use the absolute URL of the file
        lines = fh.readlines()
        #skip the vcf header that oncotator wrote out into the MAF file
        for line in lines[101:]:
            field = re.split(r'\t', line)
            if int(field[302]) == 2:
                tumor_samp.append(field[16])
                fw.write(str(line))
            elif int(field[302]) == 0:
                norm_match.append(field[16])
        fh.close()
    fw.close()

#Change the sample barcode names for the final maf file
#Replace UNKNOWN barcode columns with list of column names made above
#dat = pd.read_table(newdir + "gists_somatic_oncotate.txt")
#replace the tumor and normal barcodes
snpdat=  pd.read_csv(newdir + "gistsnp_somatic_oncotate_unlabel.txt", header=0, sep="\t", dtype=object)
snpdat.Tumor_Sample_Barcode=tumor_samp
snpdat.Matched_Norm_Sample_Barcode=norm_match
snpdat.to_csv(newdir + "gistsnp_somatic_oncotate.txt", sep="\t", dtype=object, index=False)
#only use snps where read depth was equal or above 20
#   snpdat=snpdat[snpdat['read_depth'] > 19]
#   snpdat.to_csv(newdir + "gistsnp_somatic_oncotate_q20.maf.txt", sep="\t")

#column 80 is t_ref_count and column 199=allelic depth
snpdat_trunc = snpdat[snpdat.columns[range(42)+[79,199]]]
snpdat_trunc.columns.values[[42,43]]=['T_AC','T_DP']
#snpdat_trunc = snpdat.ix[:,'Hugo_Symbol':'Protein_Change']
cond = snpdat_trunc['Tumor_Sample_Barcode'].str.contains('^RK|^SUR')

snpdatkit = snpdat_trunc[~cond]
snpdatrk = snpdat_trunc[cond]

snpdatkit.to_csv(newdir + "gistsnp_somatic_oncotate_KIT.maf.txt", sep="\t", dtype=object, index=False)
snpdatrk.to_csv(newdir + "gistsnp_somatic_oncotate_RK.maf.txt", sep="\t", dtype=object, index=False)


for subdir, dirs, files in os.walk(indeldir):
#for files in os.listdir(indeldir):
    fw = open(newdir + "gistsindel_somatic_oncotate_unlabel.txt", 'w')
    fh = open(indeldir + files[1], 'r')
    lines = fh.readlines()
    fw.write(str(lines[8]))#write column names/header into the master file
    fh.close()
    tumor_samp=[]
    norm_samp=[]
    nlstr = re.compile('.*Normal.*|.*normal.*|.*NL.*|.*Norm.*')
    for file in files:
        fh = open(indeldir + file, 'r')  # use the absolute URL of the file
        lines = fh.readlines()
        for line in lines[9:]:
            field = re.split(r'\t', line)
            if field[303]=="True\n":
                if nlstr.match(field[16]):
                    norm_samp.append(field[16])
                else:
                    tumor_samp.append(field[16])
                    fw.write(str(line))
        print file
        fh.close()
    fw.close()

#Change the sample barcode names for the final maf file
#Replace UNKNOWN barcode columns with list of column names made above
#dat = pd.read_table(newdir + "gists_somatic_oncotate.txt")
#replace the tumor and normal barcodes
indeldat=  pd.read_csv(newdir + "gistsindel_somatic_oncotate_unlabel.txt", header=0, sep="\t", dtype=object)
indeldat.Tumor_Sample_Barcode=tumor_samp
indeldat.Matched_Norm_Sample_Barcode=norm_samp
indeldat.to_csv(newdir + "gistsindel_somatic_oncotate.txt", sep="\t", dtype=object, index=False)

#only use snps where read depth was equal or above 20
#indeldat=indeldat[indeldat['T_DP'] > 19]
#columns 198, 200 have AC (allelic count) and Depth
indel_trunc=indeldat[indeldat.columns[range(42)+[198,199]]]
#indel_trunc = indeldat.ix[:,'Hugo_Symbol':'Protein_Change']
cond = indel_trunc['Tumor_Sample_Barcode'].str.contains('^RK|^SUR')

indelkit = indel_trunc[~cond]
indelrk = indel_trunc[cond]



#Combining the snps and indels into a single mafs
#Change snpdatkit alt count and read depth to match indels

mutdatkit = pd.concat([snpdatkit, indelkit])
mutdatkit.to_csv(newdir + "combo_muts_KIT.maf.txt", sep="\t", dtype=object, index=False)

mutdatrk = pd.concat([snpdatrk, indelrk])
mutdatrk.to_csv(newdir + "combo_muts_RK.maf.txt", sep="\t", dtype=object, index=False)

mutdat = pd.concat([snpdat_trunc, indel_trunc])
mutdat.to_csv(newdir + "combo_muts.maf.txt", sep="\t", dtype=object, index=False)
