#! /usr/bin/env python
'''
This script calculates fractions of SNPs with iHS values above 2.0 over
genomic windows of specified size.

#Example input:

#CHROM  POS	 iHS
chr1	14548	-3.32086
chr1	14670	-2.52
chr1	19796	0.977669
chr1	19798	3.604374
chr1	29412	-0.308192
chr1	29813	2.231736
chr1	29847	0.6594
chr1	29873	-2.03918
chr1	30050	-0.113216
chr1	30097	2.0193944
chr1	30135	-0.161264
chr1	30259	0.13628
chr1	30365	-0.357767
chr1	30370	0.953858
chr1	30664	2.0124902
chr1	30723	-0.255984
chr1	30856	3.355832
chr1	30903	-3.196446
chr1	31052	2.590459
chr1	31409	-0.497963
chr1	31414	0.611446
chr1	31424	-0.700634
chr1	31758	2.262846
chr1	31841	-0.50899
chr1	31849	5.392066
chr1	31860	-0.383864
chr1	31864	6.39043
chr1	32008	0.00886538
chr1	32158	-3.451976
chr1	32360	0.194424
chr1	32439	-0.995733


#Example output:

#CHROM	POS	nSNPs	iHS
chr1	14609.0	2	1.0
chr1	19797.0	2	0.0
chr1	29642.5	4	0.5
chr1	30476.5	10	0.4
chr1	31458.0	9	0.444444444444
chr1	32223.5	4	0.25


#command:

$ python calculate_iHSproportion.py \
    -i iHS.txt \
    -o iHS.window.txt \
    -w 1000 \
    -t 2

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls  # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the input file',
    type=str,
    required=True)
parser.add_argument(
    '-o', '--output',
    help='name of the output file',
    type=str,
    required=True)
parser.add_argument(
    '-w',
    '--window',
    help='sliding window size',
    type=int,
    required=True)
parser.add_argument(
    '-t',
    '--threshold',
    help='iHS threshold to calculate propotion for',
    type=int,
    required=True)
args = parser.parse_args()

############################# functions #############################


def proportionWindow(values, threshold):
    ''' calculates proportion of a values larger than threshold'''
    largerThan = []
    for i in values:
        if abs(i) >= threshold:
            largerThan.append(i)
    windowSize = len(values)
    proportion = len(largerThan) / float(windowSize)
    return [windowSize, proportion]


############################# program #############################

print('Opening the file...')

windSize = args.window
windPosEnd = windSize
counter = 0

with open(args.input) as datafile:
    header_line = datafile.readline()

    # make output header
    header_words = header_line.split()
    chrPos = header_words[0:2]
    chrPosP = '\t'.join(str(s) for s in chrPos)
    outputFile = open(args.output, 'w')
    outputFile.write("%s\tnSNPs\t%s\n" % (chrPosP, header_words[2]))

    print('Processing the data  ...')

    Vwindow = []
    ChrPrevious = ''
    posS = ''
    posE = ''
    for line in datafile:
        words = line.split()
        Chr = words[0]
        pos = int(words[1])
        indVal = float(words[2])

        # to store the values of a previous line
        if not ChrPrevious:
            ChrPrevious = Chr
        if not posS:
            posS = pos
        if not posE:
            posE = pos

        # if window size is reached output the results
        if Chr != ChrPrevious:  # if end of a chromosome
            meanValWindow = proportionWindow(Vwindow, args.threshold)
            meanValWindowP = '\t'.join(str(s) for s in meanValWindow)
            calls.processWindow(ChrPrevious, posS, posE,
                                meanValWindowP, outputFile)
            windPosEnd = windSize
            Vwindow = []
            posS = pos
        elif pos > windPosEnd:  # if end of a window
            if Vwindow:
                meanValWindow = proportionWindow(Vwindow, args.threshold)
                meanValWindowP = '\t'.join(str(s) for s in meanValWindow)
                calls.processWindow(Chr, posS, posE,
                                    meanValWindowP, outputFile)
            windPosEnd = windPosEnd + windSize
            Vwindow = []
            posS = pos
            while pos > windPosEnd:  # gap is larger than window size
                windPosEnd = windPosEnd + windSize

        ChrPrevious = Chr
        posE = pos

        # append values
        Vwindow.append(indVal)

        # track progress
        counter += 1
        if counter % 1000000 == 0:
            print str(counter), "lines processed"

# process the last window
meanValWindow = proportionWindow(Vwindow, args.threshold)
meanValWindowP = '\t'.join(str(s) for s in meanValWindow)
calls.processWindow(Chr, posS, pos, meanValWindowP, outputFile)

datafile.close()
outputFile.close()
print('Done!')
