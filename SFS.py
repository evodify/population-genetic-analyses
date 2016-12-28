#! /usr/bin/env python
'''
This script reconstructs the site frequency spectrum for a given set of samples.

input.file:

#CHROM xPOS sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9 sample10 sample11 sample12 sample13 sample14 sample15 sample16 sample17 sample18 sample19 sample20 sample21 sample22
scaffold_1 585 C C C A C N N C C C C C C C C C C C C C C C
scaffold_1 586 G G G G G G G G G G G G G G G G G G G G G G
scaffold_1 587 G G C C G G G G G G G G G G G G G G G G G G
scaffold_1 589 C C C C C C C C C C C C C C C C C C C C C C
scaffold_1 591 G G G G G G G G G G G G G G G G G G G G G G
scaffold_1 593 T T T T T T A A A A A A A A A A A A A A A A
scaffold_1 594 T T T T G G G G G G G G G G G G G G G G G G
scaffold_1 595 G G G G G G G G G G G G G G G G G G G G G G
scaffold_1 596 C C C C C C C C C C C C C C C C C C C C C C
scaffold_1 597 T T T T T T T T T T T T T T T T T T T T T T

ancestor.file:

#CHROM  xPOS    Ancestor
scaffold_1  585 A
scaffold_1  586 G
scaffold_1  587 G
scaffold_1  589 C
scaffold_1  591 G
scaffold_1  593 T
scaffold_1  594 T
scaffold_1  595 G
scaffold_1  596 T
scaffold_1  597 A

Output:

Derived-allele-freq Number-of-sites
0   4
2   1
16  1
18  1
22  2

The histograms of SFS is also produced.

# command

$ python SFS.py -i input.file -o output.file -a ancestor.file

contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu
'''

############################# modules #############################

import calls # my custom module
import collections
import matplotlib
matplotlib.use('Agg') # to avoid RuntimeError('Invalid DISPLAY variable'). Must be before importing matplotlib.pyplot!
import matplotlib.pyplot as plt
import numpy as np

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-a', '--ancestor', help = 'name of the file with ancestral sequence to polarize alleles', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples to process (optional)', type=str, required=False)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

ref = open(args.ancestor, 'r')
ref_header = ref.readline()
words2 = ref.readline().split()
ref_chr_pos = words2[0:2]
ref_ch = int(ref_chr_pos[0].split('_')[1])
ref_pos = int(ref_chr_pos[1])
ancest = words2[2]

output = open(args.output, 'w')

counter = 0
freq = []
with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()

  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)

  for line in datafile:
    words = line.split()
    ch = int(words[0].split('_')[1])
    pos = int(words[1])

    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"
      
    # select samples
    genotypes = calls.selectSamples(sampCol, words)

    # skip sites with missing data
    if 'N' in genotypes:
      continue
    numAl = collections.Counter(genotypes)
    numAlM = numAl.most_common()

    while ch > ref_ch or (ch == ref_ch and pos > ref_pos):
      words2 = ref.readline().split()
      if words2 == []:
        ancest = 'N'
        break
      else:
        ref_chr_pos = words2[0:2]
        ref_ch = int(ref_chr_pos[0].split('_')[1])
        ref_pos = int(ref_chr_pos[1])
        ancest = words2[2]

    if ancest == 'N':  # skip unpolarized sites
      continue
    elif len(numAl) == 1:  # if fixed
      al1 = numAlM[0][0]
      if al1 == ancest:
        freq.append(0)
      else:
        x1 = numAlM[0][1]
        freq.append(x1)
    elif len(numAl) == 2:  # if biallelic
      al1 = numAlM[0][0]
      al2 = numAlM[1][0]
      if al2 == ancest:
        x1 = numAlM[0][1]
        freq.append(x1)
      elif al1 == ancest:
        x2 = numAlM[1][1]
        freq.append(x2)

  SFSnum = collections.Counter(freq)
  SFSnumS = collections.OrderedDict(sorted(SFSnum.items()))
  output.write("Derived-allele-freq\tNumber-of-sites\n")
  for key, value in SFSnumS.items():
    output.write("%s\t%s\n" % (key, value))

# Plot frequency
bins = np.arange(len(sampCol)+1)[1:]-0.5
heightF = 6
# increase figure width if the number of samples is too large
if len(sampCol) > 50:
  widthF = float(len(sampCol))/5.0
else:
  widthF = 8
plt.figure(figsize=(widthF, heightF))
plt.hist(freq, color="grey", bins=bins)
plt.xlim(0.5, len(sampCol)-0.5)
plt.xticks(range(1,len(sampCol), 1))
plt.ylabel("Number of sites")
plt.xlabel("Derived allele frequency")
plt.title(args.output, size = 18)
plt.tight_layout()
plt.savefig(args.output+".png", dpi=90)
plt.close()

datafile.close()
output.close()
ref.close()
print('Done!')
