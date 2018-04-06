#! /usr/bin/env python

'''
This script calculates total mean of a table with numbers.

#Example input:

CHROM   POS sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   15  3   0   5   18  9   2   13
chr_1   2   13  3   0   4   16  6   4   9
chr_1   3   14  1   5   3   3   6   5   2
chr_1   4   14  1   5   3   3   7   5   2
chr_1   5   15  4   8   9   20  11  8   9
chr_2   1   9   2   13  2   1   1   2   10
chr_2   2   9   2   12  2   1   1   3   10
chr_2   3   10  3   12  2   14  1   3   7
chr_2   4   11  3   11  2   13  1   3   7


#Example output:

input.tab


#command:

$ python calculate_TotalMean.py -i input.tab -o output.tab -s "sample1,sample2,sample3,sample8"

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import calls # my custom module
import numpy as np
import warnings

############################# options #############################

parser = calls.CommandLineParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()


# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# functions #############################


############################# program #############################

print('Opening the file...')

counter = 0

with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  # count number of sample
  nSample = len(sampleNames)
############################## perform counting ####################

  print('Calculating ...')
  sumT = 0
  sumN = 0

  for line in datafile:
    words = line.split()

    # select samples
    sample_num = calls.selectSamples(sampCol, words)

    # sum up
    for i in sample_num:
      if i.isdigit():
        sumT += float(i)
        sumN += 1
      elif i == "NA":
        continue
      else:
        warnings.warn("%s is not numeric at the line %s" % (i, counter+1))

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

# make output header
outputFile = open(args.output, 'w')
heteroT = round(sumT/sumN, 2)
outputFile.write("%s\t%s\n" % (args.input, heteroT))

datafile.close()
outputFile.close()
print('Done!')

