oncohaplodir = "C:\Users\Owner\Box Sync\oncotator\oncotator_haplo\\"
newdir = "C:\Users\Owner\Box Sync\haplotypecaller\haplofiltered\\"


def onco_haplo(file):
    from

for subdir, dirs, files in os.walk(oncohaplodir):
    #input the directory you want to walk through full of the oncotate files
    fw = open(newdir + "gistsnp_somatic_oncotate_unlabel.txt", 'w')
    #fw = open(newdir + file.replace("_merged_sort_mutect_filt_annotated", "_somatic_oncotate"), 'w')
    fh = open(snpdir + files[1], 'r')
    lines = fh.readlines()
    fw.write(str(lines[100]))#write column names/header into the master file
    fh.close()
    norm_match=[]
    tumor_samp=[]
    for file in files:
        fh = open(snpdir + file, 'r')  # use the absolute URL of the file
        lines = fh.readlines()
        for line in lines[101:]:
            field = re.split(r'\t', line)
            if int(field[302]) == 2:
                tumor_samp.append(field[16])
                fw.write(str(line))
            elif int(field[302]) == 0:
                norm_match.append(field[16])
        fh.close()
    fw.close()