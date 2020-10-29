#! /usr/bin/env python
'''
Calculates average values for sliding windows centred at each SNP.

#Example input:

CHROM	POS	stats1	stats2
chr1	14670	0.69	1.7
chr1	19796	0.83	1.8
chr1	29099	0.90	1.9
chr1	30259	0.80	1.8
chr1	31671	0.44	1.4
chr1	41386	0.56	1.6
chr1	50901	0.57	1.6
chr2	947	0.53	1.5
chr2	1017	0.49	1.5
chr2	51111	0.36	1.4
chr2	75365	0.67	1.7
chr2	80686	0.50	1.5
chr2	173083	0.75	1.8
chr2	173098	0.77	1.8

#Example output:

CHROM	POS	stats1	stats2
chr1	14670	0.732	1.72
chr1	19796	0.703333333333	1.7
chr1	29099	0.684285714286	1.68571428571
chr1	30259	0.684285714286	1.68571428571
chr1	31671	0.684285714286	1.68571428571
chr1	41386	0.683333333333	1.68333333333
chr1	50901	0.654	1.66
chr2	947	0.51	1.5
chr2	1017	0.51	1.5
chr2	51111	0.515	1.55
chr2	75365	0.51	1.53333333333
chr2	80686	0.585	1.6
chr2	173083	0.76	1.8
chr2	173098	0.76	1.8

#command:

$ python calculate_WindowAverage_perSNP.py \
    -i input.tab \
    -o output.tab \
    -w 50000

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls  # my custom module
import bisect

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i', '--input', help='name of the input file', type=str, required=True)
parser.add_argument(
    '-o', '--output', help='name of the output file', type=str, required=True)
parser.add_argument(
    '-w', '--window', help='sliding window size', type=int, required=True)
parser.add_argument(
    '-m', '--min', help='minimum number of SNPs to keep window',
    type=int, required=False, default=1)
args = parser.parse_args()

############################# functions #############################

def average(lst, minNumberSNPs):
    if len(lst) >= minNumberSNPs:
        avg = sum(lst)/len(lst)
    else:
        avg = 'NA'
    return avg

def slideWindow():
    for i in posList:
        # if start of chromosome
        if i <= windSize:
            startIndex = 0
        else:
            windStart = i-windSize
            startIndex = bisect.bisect(posList, windStart)
        windowEnd = i+windSize
        endIndex = bisect.bisect(posList, windowEnd)

        outputFile.write("%s\t%s" % (ChrPrevious, i))
        for n in statName:
            windowAverage = average(statsDic[n][startIndex:endIndex], minSNPs)
            outputFile.write("\t%s" % windowAverage)
        outputFile.write("\n")

############################# program #############################

print('Opening the file...')

windSize = args.window/2
counter = 0
minSNPs = args.min 

with open(args.input) as datafile:
    header_line = datafile.readline()
    statName = header_line.split()[2:]

    # make output header
    outputFile = open(args.output, 'w')
    outputFile.write(header_line)

    print('Processing the data  ...')

    ChrPrevious = ''
    posStart = ''
    posEnd = ''
    for line in datafile:
        words = line.split()
        Chr = words[0]
        pos = int(float(words[1]))
        stats = [float(i) for i in words[2:]]
        
        # read in one chromosome
        if ChrPrevious == '':
            posList = []
            statsDic = {key: [] for key in statName}
            posList.append(pos)
            for n in statName:
                statsDic[n].append(stats[statName.index(n)])
        elif Chr == ChrPrevious:
            posList.append(pos)
            for n in statName:
                statsDic[n].append(stats[statName.index(n)])
        else:
            slideWindow()
            posList = []
            statsDic = {key: [] for key in statName}
            posList.append(pos)
            for n in statName:
                statsDic[n].append(stats[statName.index(n)])
           
        ChrPrevious = Chr
        # track progress
        counter += 1
        if counter % 1000000 == 0:
            print str(counter), "lines processed"

    # process the last chromosome
    slideWindow()

datafile.close()
outputFile.close()
print('Done!')
