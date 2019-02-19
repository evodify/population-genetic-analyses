#! /usr/bin/env python
'''
This script calculates average values per sliding window.

#Example input:

CHROM  POS  sample1  sample2  sample3  sample4  sample5
chr1    2923    0       16      13      24      27
chr1    4696    1       3       5       13      6
chr1    6240    5       10      5       15      19
chr1    6244    5       10      5       16      20
chr1    6527    9       20      12      20      36
chr1    6544    NA       21      16      20      36
chr1    6665    5       17      12      15      32
chr1    6676    5       22      14      18      31
chr1    6677    5       22      14      18      31
chr1    8017    14      19      9       20      33
chr1    8374    12      5       16      13      24
chr1    8618    7       13      10      25      21
chr1    8986    16      19      10      34      20
chr1    9185    15      31      18      42      44
chr1    9218    15      30      21      45      45
chr1    9374    16      28      18      45      43
chr1    9378    16      27      19      43      42
chr1    9411    18      24      NA      50      42
chr1    10743   10      17      16      34      28
chr1    11105   47      36      46      66      69
chr1    11162   14      24      32      43      55
chr1    11331   45      34      82      41      87
chr1    11368   51      41      107     57      101
chr1    13956   17      15      33      38      32
chr1    14548   5       4       10      9       8
chr1    14670   22      16      51      NA      22
chr1    14686   22      35      57      63      42
chr1    19796   54      32      43      57      49
chr1    19798   54      32      45      56      48


#Example output:

CHROM  POS  sample1  sample2  sample3  sample4  sample5
chr1	3809	1	10	9	19	17
chr1	7825	11	20	13	27	32
chr1	12714	26	25	48	44	49
chr1	19797	54	32	44	57	49


#command:

$ python calculate_AveragePerWindow.py -i input.tab -o output.tab -w 5000

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls  # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument(
    '-i', '--input', help='name of the input file', type=str, required=True)
parser.add_argument(
    '-o', '--output', help='name of the output file', type=str, required=True)
parser.add_argument(
    '-w', '--window', help='sliding window size', type=int, required=True)
args = parser.parse_args()

############################# functions #############################

def meanWindow(dictList):
    ''' calculates mean of a window'''
    for k in dictList:
        values = []
        for val in dictList[k]:
            if val != 'NA':
                values.append(float(val))
        if len(values) > 0:
            averageValue = sum(values) / len(values)
            dictList[k] = averageValue
        else:
            dictList[k] = 'NA'
    return dictList

def createNewDict(NamesList):
    ''' creates a new empty dictionary with sample names as keys'''
    newDict = {}
    for k in NamesList:
        newDict[k] = []
    return newDict

def printWindow(inputDict, orderedNames):
    ''' creates print string from a dictionary with mean values'''
    newList = []
    for n in orderedNames:
        newList.append(inputDict[n])
    newListP = '\t'.join(str(el) for el in newList)
    return newListP

############################# program #############################

print('Opening the file...')

windSize = args.window
windPosEnd = windSize
counter = 0

with open(args.input) as datafile:
    header_line = datafile.readline()

    # make output header
    outputFile = open(args.output, 'w')
    outputFile.write(header_line)

    # make samples dict
    header_words = header_line.split()
    sampleNames = header_words[2:]
    windowDict = createNewDict(sampleNames)

    print('Processing the data  ...')

    ChrPrevious = ''
    posS = ''
    posE = ''
    for line in datafile:
        words = line.split()
        Chr = words[0]
        pos = int(float(words[1]))
        indVal = words[2:]

        # to store the values of a previous line
        if not ChrPrevious:
            ChrPrevious = Chr
        if not posS:
            posS = windPosEnd - windSize
        if not posE:
            posE = windPosEnd

        # if window size is reached output the results
        if Chr != ChrPrevious:  # if end of a chromosome
            meanValWindow = meanWindow(windowDict)
            meanValWindowP = printWindow(meanValWindow, sampleNames)
            calls.processWindow(ChrPrevious, posS, posE,
                                meanValWindowP, outputFile)
            windPosEnd = windSize
            windowDict = createNewDict(sampleNames)
            posS = windPosEnd - windSize
        elif pos > windPosEnd:  # if end of a window
            meanValWindow = meanWindow(windowDict)
            meanValWindowP = printWindow(meanValWindow, sampleNames)
            calls.processWindow(Chr, posS, posE,
                                meanValWindowP, outputFile)
            windPosEnd = windPosEnd + windSize
            windowDict = createNewDict(sampleNames)
            posS = windPosEnd - windSize
            while pos > windPosEnd:  # gap is larger than window size
                windPosEnd = windPosEnd + windSize

        ChrPrevious = Chr
        posE = windPosEnd

        # append values
        for s in xrange(len(sampleNames)):
            windowDict[sampleNames[s]].append(indVal[s])

        # track progress
        counter += 1
        if counter % 1000000 == 0:
            print str(counter), "lines processed"

# process the last window
meanValWindow = meanWindow(windowDict)
meanValWindowP = printWindow(meanValWindow, sampleNames)
calls.processWindow(Chr, posS, windPosEnd,
                    meanValWindowP, outputFile)

datafile.close()
outputFile.close()
print('Done!')
