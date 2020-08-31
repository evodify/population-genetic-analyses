#! /usr/bin/env python
'''
Calculates average values for sliding windows centred at each SNP.

#Example input:

CHROM  POS  stats
chr1	14670	0.68568
chr1	19796	0.832316
chr1	29099	0.895923
chr1	30259	0.80235
chr1	31671	0.437556999117073
chr1	41386	0.555906
chr1	50901	0.566755
chr2	947	0.528296
chr2	1017	0.494082466499707
chr2	51111	0.360812988529893
chr2	75365	0.665836
chr2	80686	0.500991
chr2	173083	0.753236
chr2	173098	0.767467

#Example output:

CHROM  POS  stats
chr1	14670	0.730765199823
chr1	19796	0.701621999853
chr1	29099	0.682355285588
chr1	30259	0.682355285588
chr1	31671	0.682355285588
chr1	41386	0.68180116652
chr1	50901	0.651698199823
chr2	947	0.51118923325
chr2	1017	0.51118923325
chr2	51111	0.513324494265
chr2	75365	0.50921332951
chr2	80686	0.5834135
chr2	173083	0.7603515
chr2	173098	0.7603515


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
    type=int, required=True, default=1)
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

        windowAverage = average(statsList[startIndex:endIndex], minSNPs)

        outputFile.write("%s\t%s\t%s\n" % (ChrPrevious, i, windowAverage))

############################# program #############################

print('Opening the file...')

windSize = args.window/2
counter = 0
minSNPs = args.min 

with open(args.input) as datafile:
    header_line = datafile.readline()

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
        stats = float(words[2])
        
        # read in one chromosome
        if ChrPrevious == '':
            posList = []
            statsList = []
            posList.append(pos)
            statsList.append(stats)
        elif Chr == ChrPrevious:
            posList.append(pos)
            statsList.append(stats)
        else:
            slideWindow()
            posList = []
            statsList = []
            posList.append(pos)
            statsList.append(stats)
           
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
