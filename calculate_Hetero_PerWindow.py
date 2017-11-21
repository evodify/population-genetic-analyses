#! /usr/bin/env python

'''
This script calculates heterozygosity with the sliding window approach.

#Example input:

CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   N   A   N   N   N   N
chr_1   2   C   Y   Y   N   C   C   N   C   N
chr_1   3   C   N   C   N   C   C   C   C   C
chr_1   4   T   T   T   N   T   T   T   T   T
chr_2   1   A   A   A   N   A   A   A   A   A
chr_2   2   C   C   C   N   C   C   C   C   C
chr_2   3   C   N   N   N   N   N   N   N   N
chr_2   4   C   C   T   C   C   C   C   C   C
chr_2   5   T   T   C   T   Y   T   Y   T   T
chr_3   1   G   G   N   N   G   N   N   N   N
chr_3   2   C   S   C   N   C   C   N   C   N
chr_3   3   N   N   N   N   N   N   N   N   N
chr_3   4   N   T   T   N   T   T   T   T   N
chr_3   5   G   -   N   N   G   G   G   C   G


#Example input2:

CHROM POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1 1 A/A A/T ./. ./. A/A ./. ./. ./. ./.
chr_1 2 C/C T/C T/C ./. C/C C/C ./. C/C ./.
chr_1 3 C/C ./. C/C ./. C/C C/C C/C C/C C/C
chr_1 4 T/T T/T T/T ./. T/T T/T T/T T/T T/T
chr_2 1 A/A A/A A/A ./. A/A A/A A/A A/A A/A
chr_2 2 C/C C/C C/C ./. C/C C/C C/C C/C C/C
chr_2 3 C/C ./. ./. ./. ./. ./. ./. ./. ./.
chr_2 4 C/C C/C T/T C/C C/C C/C C/C C/C C/C
chr_2 5 T/T T/T C/C T/T T/C T/T T/C T/T T/T
chr_3 1 G/G G/G ./. ./. G/G ./. ./. ./. ./.
chr_3 2 C/C G/C C/C ./. C/C C/C ./. C/C ./.
chr_3 3 ./. ./. ./. ./. ./. ./. ./. ./. ./.
chr_3 4 ./. T/T T/T ./. T/T T/T T/T T/T ./.
chr_3 5 G/G -/- ./. ./. G/G G/G G/G C/C G/G

#Example output:

CHROM   POS Ns
chr_1   2.5 3.0
chr_1   6.0 2.0
chr_2   3.0 2.0
chr_3   3.0 4.2

#command:

$ python calculate_Hetero_PerWindow.py -i input.tab -o output.tab -w 5 -s "sample1,sample2,sample3,sample4,sample5,sample6,sample7,sample8"

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
parser.add_argument('-w', '--window', help = 'sliding window size', type=int, required=True)
parser.add_argument('-m', '--missing', help = 'number of allowed Ns per position ', type=int, required=False)
args = parser.parse_args()


# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

if args.missing:
 allowedN = args.missing
else:
   allowedN = 0

############################# functions #############################

def meanWindow(hetero, total):
  ''' calculates mean of a window'''
  if total:
    propHetero = sum(hetero)/sum(total)
  else:  # is all sites are missing data
    propHetero = 'NA'
  return propHetero

############################# program #############################

print('Opening the file...')

windSize = args.window
windPosEnd = windSize
counter = 0

with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  # count number of sample
  nSample = len(sampleNames)

  # make output header
  outputFile = open(args.output, 'w')
  outputFile.write("CHROM\tPOS\tHeter\n")

############################## perform counting ####################

  print('Counting heterozygots ...')

  Hwindow = []
  Twindow = []
  ChrPrevious = ''
  posS = ''
  posE = ''
  for line in datafile:
    words = line.split()
    Chr = words[0]
    pos = int(words[1])

    # to store the values of a previous line
    if not ChrPrevious:
      ChrPrevious = Chr
    if not posS:
      posS = pos
    if not posE:
      posE = pos

    # select samples
    sample_charaters = calls.selectSamples(sampCol, words)

    # check if one- or two-character code
    if any(["/" in gt for gt in sample_charaters]):
      sample_charaters = calls.twoToOne(sample_charaters)

    # if window size is reached output the results
    if Chr > ChrPrevious:  # if end of a chromosome
      try:
        HeterWindow = round(meanWindow(Hwindow, Twindow), 4)
      except Exception:
        HeterWindow = "NA"
      calls.processWindow(ChrPrevious, posS, posE, HeterWindow, outputFile)
      windPosEnd = windSize
      Hwindow = []
      Twindow = []
      posS = pos
    elif pos > windPosEnd:  # if end of a window
      try:
        HeterWindow = round(meanWindow(Hwindow, Twindow), 4)
      except Exception:
        HeterWindow = "NA"
      calls.processWindow(Chr, posS, posE, HeterWindow, outputFile)
      windPosEnd = windPosEnd+windSize
      Hwindow = []
      Twindow = []
      posS = pos
      while pos > windPosEnd:  # if the gap in positions is larger than window size
        windPosEnd = windPosEnd+windSize

    ChrPrevious = Chr
    posE = pos

    # count hetero
    Nmising = calls.countPerPosition(sample_charaters, 'N')
    if Nmising < allowedN: # skip if too many Ns
      nHerer = calls.countHeteroPerPosition(sample_charaters)
      nTotal = float(nSample - Nmising)
      Hwindow.append(float(nHerer))
      Twindow.append(float(nTotal))

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

# process the last window
try:
  HeterWindow = round(meanWindow(Hwindow, Twindow), 4)
except Exception:
  HeterWindow = "NA"
calls.processWindow(Chr, posS, pos, HeterWindow, outputFile)

datafile.close()
outputFile.close()
print('Done!')

