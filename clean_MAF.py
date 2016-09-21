#!/usr/bin/env python
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Output MAF file from combining Oncotator output, which had combined tumor normals per line
Usage:
import clean_MAF
snpdat=clean_MAF.snps2MAF()
indeldat=clean_MAF.indel2MAF()
clean_MAF.trunc2comb(snpdat, indeldat)
#Less modules Indel and SNP mutect output are very completely different
#because I'm cursed. Something fixable, when I'm smarter.
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
import os
import re
import pandas as pd

newdir = "C:\Users\Owner\Box Sync\mutsig\inputs\\"
snpdir = "C:\\Users\\Owner\\Box Sync\\oncotator\\oncotator_snps\\"
indeldir = "C:\Users\Owner\Box Sync\oncotator\oncotator_indel\\"


def snps2MAF():
    with open(snpdir + os.listdir(snpdir)[1], 'r') as fh:
        lines = fh.readlines()
        headerline = re.split(r'\t', lines[100])
    norm_name = []
    tumor_samp = []
    tumor_name = []
    for file in os.listdir(snpdir):
        fh = open(snpdir + file, 'r')  # use the absolute URL of the file
        lines = fh.readlines()
        for line in lines[101:]:  # skip the vcf header that oncotator wrote out into the MAF file
            field = re.split(r'\t', line)
            if int(field[302]) == 2:
                tumor_samp.append(field)
                tumor_name.append(field[16])
            elif int(field[302]) == 0:
                norm_name.append(field[16])
        fh.close()
    snpdat = pd.DataFrame(tumor_samp, columns=headerline)
    # Change the sample barcode names for the final maf file
    # Replace UNKNOWN barcode columns with list of column names made above
    snpdat.Tumor_Sample_Barcode = tumor_name
    snpdat.Matched_Norm_Sample_Barcode = norm_name
    snpdat.to_csv(newdir + "gistsnp_somatic_oncotate.txt", sep="\t", dtype=object, index=False)
    return snpdat

# only use snps where read depth was equal or above 20
# snpdat=snpdat[snpdat['read_depth'] > 19]
# snpdat.to_csv(newdir + "gistsnp_somatic_oncotate_q20.maf.txt", sep="\t")


def indel2MAF():
    with open(indeldir + os.listdir(indeldir)[1], 'r') as fh:
        lines = fh.readlines()
        headerline = re.split(r'\t', lines[8])
    tumor_samp = []
    tumor_name = []
    norm_name = []
    nlstr = re.compile('.*Normal.*|.*normal.*|.*NL.*|.*Norm.*')
    for file in os.listdir(indeldir):
        fh = open(indeldir + file, 'r')  # use the absolute URL of the file
        lines = fh.readlines()
        for line in lines[9:]:
            field = re.split(r'\t', line)
            if field[303] == "True\n":
                if nlstr.match(field[16]):
                    norm_name.append(field[16])
                else:
                    tumor_name.append(field[16])
                    tumor_samp.append(field)
        fh.close()

    # Change the sample barcode names for the final maf file
    # Replace UNKNOWN barcode columns with list of column names made above
    indeldat = pd.DataFrame(tumor_samp, columns=headerline)
    indeldat.Tumor_Sample_Barcode = tumor_name
    indeldat.Matched_Norm_Sample_Barcode = norm_name
    indeldat.to_csv(newdir + "gistsindel_somatic_oncotate.txt", sep="\t", dtype=object, index=False)
    return indeldat
# only use snps where read depth was equal or above 20
# indeldat=indeldat[indeldat['T_DP'] > 19]
# columns 198, 200 have AC (allelic count) and Depth


def trunc2comb(snpdat, indeldat):
    """
    Combining the snps and indels into a single mafs
    :param snpdat:
    :param indeldat:
    :return:nothing output is files under newdir for MutSig analysis
    """

    # Change snpdatkit alt count and read depth to match indels
    # column 80 is t_ref_count and column 199=allelic depth
    snpdat_trunc = snpdat[snpdat.columns[range(42) + [300, 199]]]
    snpdat_trunc.columns.values[[42, 43]] = ['T_DP', 'T_AC']

    condsnp = snpdat_trunc['Tumor_Sample_Barcode'].str.contains('^RK|^SUR')
    snpdatkit = snpdat_trunc[~condsnp]
    snpdatrk = snpdat_trunc[condsnp]

    # columns 198, 200 have AC (allelic count) and Depth
    indel_trunc = indeldat[indeldat.columns[range(42)+[199, 198]]]
    condindel = indel_trunc['Tumor_Sample_Barcode'].str.contains('^RK|^SUR')
    indelkit = indel_trunc[~condindel]
    indelrk = indel_trunc[condindel]
    mutdatkit = pd.concat([snpdatkit, indelkit])
    mutdatkit.to_csv(newdir + "combo_muts_KIT.maf.txt", sep="\t", dtype=object, index=False)

    mutdatrk = pd.concat([snpdatrk, indelrk])
    mutdatrk.to_csv(newdir + "combo_muts_RK.maf.txt", sep="\t", dtype=object, index=False)

    mutdat = pd.concat([snpdat_trunc, indel_trunc])
    mutdat.to_csv(newdir + "combo_muts.maf.txt", sep="\t", dtype=object, index=False)
    print "all done"

# only use snps where read depth was equal or above 20 with commented code below
# indeldat=indeldat[indeldat['T_DP'] > 19] # columns 198, 200 have AC (allelic count) and Depth




