#! /usr/bin/env python
'''
This script calculates change in derived allele frequency (deltaDAF) between
two populations in a polarized tab-delimited file.
The selected population should be the first one in the samples field.

# input.file:
CHROM	POS	REF	ALT	sample1	sample2	sample3	sample4	sample5	sample6
chr6	215	T	G	1|0	1|0	0|2	1|0	0|0	0|0
chr6	245	A	C	0|0	0|0	1|1	1|1	1|1	1|1
chr6	261	C	A	0|0	0|0	0|0	0|0	0|0	0|0
chr6	280	C	A	0|0	1|0	0|0	0|0	0|0	0|0
chr6	321	G	C	0|0	0|0	0|0	0|0	0|0	0|0

# output:
CHROM	POS	ANC	DER	pop1Freq	pop2Freq	deltaDAF
chr6	245	.	A	0.33333	1.0	-0.66667
chr6	280	.	C	0.16667	0.0	0.16667

# command:

$ python DAF.py -i input.file -o output.file \
    -s "pop1[sample1,sample2,sample3];pop2[sample4,sample5,sample6]"

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls  # my custom module
import re # to split input
from collections import Counter # for counting
############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i', '--input',
    help='name of the input file',
    type=str, required=True)
parser.add_argument(
    '-o', '--output',
    help='name of the output file',
    type=str,
    required=True)
parser.add_argument(
    '-s',
    '--samples',
    help='Sample names with first population being testes as selected one:\
        -s "pop1[sample1,...];pop2[sample4,...]"',
    type=str,
    required=True)
args = parser.parse_args()

############################# program #############################

with open(args.input) as datafile:
    # process header
    header = datafile.readline()
    while header.startswith("##"):
        header = datafile.readline()
    header_words = header.split()

    # check if all sample names are present in the header
    popNames = args.samples
    selPop = popNames.split("[")[0] # define the selected population
    Psamples = []
    popSamples = {}
    populations = popNames.strip("\"").split(";")
    if len(populations) != 2:
        raise IOError(
                    'Exactly two populations are required.\
                     You provided "%s" populations' % len(populations))
    for i in populations:
        pName = i.split("[")[0]
        pSample = re.split("\[|\]|", i)[1]
        Psamples = (pSample.split(","))
        for s in Psamples:
            if s not in header_words:
                raise IOError(
                    'Sample name "%s" is not found in the header' % (s))
        popSamples[pName] = Psamples
    
    # create the output file
    output = open(args.output, 'w')
    print('Creating the output file...')
    selPopP = str(selPop)+"Freq"
    pNameP = str(pName)+"Freq"
    output.write("CHROM\tPOS\tANC\tDER\t%s\t%s\tdeltaDAF\n" % (selPopP, pNameP))


    # index samples
    popCol = {}
    for pop in popSamples:
        popCol[pop] = calls.indexSamples(popSamples[pop], header_words)

    for line in datafile:
        words = line.split()
        chr = words[0]
        pos = int(words[1])
        anc = words[2]
        der = words[3]

        # calculate frequencies of derived allele
        popFreq = {}
        for pop in popCol:
            gt = calls.selectSamples(popCol[pop], words)
            gtSplit = []
            for g in gt:
                gtSplit.append(g.split("|"))
            gtSplitF = calls.flattenList(gtSplit)
            if len(set(gtSplitF)) <= 2:
                gtCount = Counter(gtSplitF)
                popFreq[pop] = float(gtCount['1'])/float(len(gtSplitF))
            else:
                popFreq[pop] = "NA"
        
        # calculate delta DAF
        if (popFreq[selPop] != 'NA' and popFreq[pName] != 'NA') and \
           (popFreq[selPop] != 0.0 or popFreq[pName] != 0.0):
            deltaDAF = popFreq[selPop] - popFreq[pName]
            deltaDAFP= round(deltaDAF, 5)
            Pop1FreqP= round(popFreq[selPop], 5)
            Pop2FreqP= round(popFreq[pName], 5)

            output.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr, pos, anc, der, 
                         Pop1FreqP, Pop2FreqP, deltaDAFP))       

datafile.close()
output.close()
print('Done!')
