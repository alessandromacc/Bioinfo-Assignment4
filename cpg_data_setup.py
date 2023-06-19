import random
from matplotlib import pyplot

def seqExtraction(annotationfile, sequencefile, newname):
    af = open(annotationfile, 'r')
    sf = open(sequencefile, 'r')
    nf = open(newname, 'w')
    sequence = ''.join(i.strip('\n').upper() for i in sf.readlines()[1:])
    for line in af.readlines():
        start = int(line.split('\t')[1])
        end = int(line.split('\t')[2])
        line_to_write = sequence[start:end]+'\n'
        nf.write(line_to_write)
    af.close()
    sf.close()
    nf.close()

#seqExtraction('chr22_annot.txt', 'chr22.fa', 'CpG.txt')

def correctLocations(af, nf):
    annot = open(af, 'r')
    newf = open(nf, 'a')
    lines = annot.readlines()
    ref = int(lines[0].split('\t')[1])
    for i in lines:
        line = i.split('\t')
        line[1] = str(int(line[1]) - ref)
        line[2] = str(int(line[2]) - ref)
        newf.write('\t'.join(line))
    annot.close()
    newf.close()

#correctLocations('chr22_annot.txt', 'chr22_annot_corrected.txt')

def generateRandomSequences(fasta, cpgfile, nf):
    sf = open(fasta, 'r')
    cpgislnds = open(cpgfile, 'r')
    outside = open(nf, 'a')

    sequence = ''.join(i.strip('\n').upper() for i in sf.readlines()[1:])
    gsize = len(sequence)
    for line in cpgislnds.readlines():
        line = line.strip('\n')
        lsize = len(line)
        rsp = random.randint(0, gsize-1-lsize)
        if 'N' not in sequence[rsp:rsp+lsize]:
            outside.write(sequence[rsp:rsp+lsize])
            outside.write('\n')
        else:
            while 'N' in sequence[rsp:rsp+lsize]:
                rsp = random.randint(0, gsize-1-lsize)
                if 'N' not in sequence[rsp:rsp+lsize]:
                    outside.write(sequence[rsp:rsp+lsize])
                    outside.write('\n')
    sf.close()
    cpgislnds.close()
    outside.close()

#generateRandomSequences('chr22.fa', 'CpG.txt', 'outside.txt')

def plot_annotated_cpg(af):
    annot = open(af, 'r')
    sites = []
    for line in annot.readlines():
        line = line[:-1].split('\t')
        sites.append(line[1])
    annot.close()
    pyplot.plot()
    for i in sites:
        pyplot.axvline(i, color = 'red', linewidth = 0.4)
    pyplot.show()

plot_annotated_cpg('chr22_annot.txt')