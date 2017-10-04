#! /usr/bin/env python

'''
This script calculates average values per sliding window.

#Example input:

#CHROM  POS Value
scaffold_1  1    8
scaffold_1  2   8
scaffold_1  3   2
scaffold_1  4   10
scaffold_1  5   10
scaffold_1  6   3
scaffold_1  7   6
scaffold_1  8   4
scaffold_1  9   7
scaffold_1  10   1
scaffold_1  11   4
scaffold_1  12   3
scaffold_1  13   8
scaffold_1  14   15
scaffold_1  15   3
scaffold_1  16   9
scaffold_1  17   2
scaffold_1  18   15
scaffold_1  19   1
scaffold_1  20   9


#Example output:

#CHROM  POS Value
scaffold_1  3.0 7.6
scaffold_1  8.0 4.2
scaffold_1  13.0    6.6
scaffold_1  18.0    7.2


#command:

$ python calculate_AveragePerWindow.py -i input.tab -o output.tab -w 5

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-w', '--window', help = 'sliding window size', type=int, required=True)
args = parser.parse_args()

############################# functions #############################

def meanWindow(values):
  ''' calculates mean of a window'''
  averageValue = sum(values)/len(values)
  return averageValue

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
    if Chr > ChrPrevious:  # if end of a chromosome
      meanValWindow = meanWindow(Vwindow)
      calls.processWindow(ChrPrevious, posS, posE, meanValWindow, outputFile)
      windPosEnd = windSize
      Vwindow = []
      posS = pos
    elif pos > windPosEnd:  # if end of a window
      meanValWindow = meanWindow(Vwindow)
      calls.processWindow(Chr, posS, posE, meanValWindow, outputFile)
      windPosEnd = windPosEnd+windSize
      Vwindow = []
      posS = pos
      while pos > windPosEnd:  # if the gap in positions is larger than window size
        windPosEnd = windPosEnd+windSize

    ChrPrevious = Chr
    posE = pos

    # append values
    Vwindow.append(indVal)

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

# process the last window
meanValWindow = meanWindow(Vwindow)
calls.processWindow(Chr, posS, pos, meanValWindow, outputFile)

datafile.close()
outputFile.close()
print('Done!')

