#!/usr/bin/env python
"""
This script estimates a mutation matrix required for LDhelmet (http://dx.doi.org/10.1371/journal.pgen.1003090).

# input:

#CHROM POS sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8 sample9 sample10 sample11 sample12 sample13 sample14 sample15 sample16 sample17 sample18 sample19 sample20 sample21
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

# ancestor
#CHROM  POS    Ancestor
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

# output:

A       T       G       C
16      31      99      62

0.00    0.54    0.00    0.46
0.20    0.00    0.22    0.28
0.00    0.00    0.99    0.01
0.00    0.00    0.00    1.00

# command:

$ python mutMatrix.py -i datafile -o outputfile -s "sample1,sample2,sample3,sample4"


contact Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

#import collections
import calls # my custom module
import collections

############################# options #############################

parser = calls.MyParser()
parser.add_argument('-i', '--input', help = 'name of the input file', type=str, required=True)
parser.add_argument('-a', '--ancestral', help = 'name of the file with ancestral alleles', type=str, required=True)
parser.add_argument('-o', '--output', help = 'name of the output file', type=str, required=True)
parser.add_argument('-s', '--samples', help = 'column names of the samples for with to calculate Ns', type=str, required=False)
parser.add_argument('-N', '--missing', help = 'number of allowed Ns', type=int, required=True)
args = parser.parse_args()

# check if samples names are given and if all sample names are present in a header
sampleNames = calls.checkSampleNames(args.samples, args.input)

############################# program #############################

counter = 0

ref = open(args.ancestral, 'r')
ref_header = ref.readline()
words2 = ref.readline().split()
ref_chr_pos = words2[0:2]
ref_ch = int(ref_chr_pos[0].split('_')[1])
ref_pos = int(ref_chr_pos[1])
ancest = words2[2]

totalA = 0
totalT = 0
totalG = 0
totalC = 0

totalAT = 0
totalAG = 0
totalAC = 0

totalTA = 0
totalTG = 0
totalTC = 0

totalCT = 0
totalCG = 0
totalCA = 0

totalGT = 0
totalGA = 0
totalGC = 0

print('Opening the file...')
with open(args.input) as datafile:
  header_line = datafile.readline()
  header_words = header_line.split()
  
  # index samples
  sampCol = calls.indexSamples(sampleNames, header_words)
    
  for line in datafile:
    # track progress
    counter += 1
    if counter % 1000000 == 0:
      print str(counter), "lines processed"

    words = line.split()
    chr_pos = words[0:2]
    ch = int(words[0].split('_')[1])
    pos = int(words[1])

     # select samples
    alleles = calls.selectSamples(sampCol, words)

    # count Ns
    valueN = calls.countPerPosition(alleles, 'N')

    if valueN <= args.missing:
      Allalleles = [i for i in alleles if i != 'N']
    else:
      continue

    # find overlap with the ancestor
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
    
    # skip unpolarized sites
    if ancest == 'N':
         continue
 
    # count alleles
    numAl = collections.Counter(Allalleles)
    numAlM = numAl.most_common()
    totalA += int(numAl['A'])
    totalT += int(numAl['T'])
    totalG += int(numAl['G'])
    totalC += int(numAl['C'])
    
    n = len(Allalleles)
     
    if len(numAlM) > 2:  # skip non-biallelic
      continue
    elif chr_pos == ref_chr_pos:
      valueAnces = int(numAl[ancest])
      valueDer = int(n-valueAnces)
      AlelleDer = list(set(list(numAl))-set(ancest))
      if ancest == 'A' and AlelleDer == ['T']:
        totalAT += valueDer
      elif ancest == 'A' and AlelleDer == ['G']:
        totalAG += valueDer
      elif ancest == 'A' and AlelleDer == ['C']:
        totalAC += valueDer
      elif ancest == 'T' and AlelleDer == ['A']:
        totalTA += valueDer
      elif ancest == 'T' and AlelleDer == ['G']:
        totalTG += valueDer
      elif ancest == 'T' and AlelleDer == ['C']:
        totalTC += valueDer          
      elif ancest == 'G' and AlelleDer == ['A']:
        totalGA += valueDer
      elif ancest == 'G' and AlelleDer == ['T']:
        totalGT += valueDer
      elif ancest == 'G' and AlelleDer == ['C']:
        totalGC += valueDer
      elif ancest == 'C' and AlelleDer == ['A']:
        totalCA += valueDer
      elif ancest == 'C' and AlelleDer == ['T']:
        totalCT += valueDer
      elif ancest == 'C' and AlelleDer == ['G']:
        totalCG += valueDer

sA = float(totalAT)/float(totalA) + float(totalAG)/float(totalA) + float(totalAC)/float(totalA)
sT = float(totalTA)/float(totalT) + float(totalTG)/float(totalT) + float(totalTC)/float(totalT)
sG = float(totalGT)/float(totalG) + float(totalGA)/float(totalG) + float(totalGC)/float(totalG)
sC = float(totalCT)/float(totalC) + float(totalCG)/float(totalC) + float(totalCA)/float(totalC)
M = max(sA, sT, sG, sC)

pAT = float(totalAT) / (float(totalA) * M)
pAG = float(totalAG) / (float(totalA) * M)
pAC = float(totalAC) / (float(totalA) * M)
pAA = 1 - (pAT + pAG + pAC)
pTA = float(totalTA) / (float(totalT) * M)
pTG = float(totalTG) / (float(totalT) * M)
pTC = float(totalTC) / (float(totalT) * M)
pTT = 1 - (pAT + pAG + pAC)
pGA = float(totalGA) / (float(totalG) * M)
pGT = float(totalGT) / (float(totalG) * M)
pGC = float(totalGC) / (float(totalG) * M)
pGG = 1 - (pGT + pGA + pGC)
pCA = float(totalCA) / (float(totalC) * M)
pCG = float(totalCG) / (float(totalC) * M)
pCT = float(totalCT) / (float(totalC) * M)
pCC = 1 - (pCT + pCG + pCA)

output = open(args.output, 'w')
output.write("%.2f\t%.2f\t%.2f\t%.2f\n" % (pAA, pAT, pAG, pAC))
output.write("%.2f\t%.2f\t%.2f\t%.2f\n" % (pTA, pTT, pTG, pTC))
output.write("%.2f\t%.2f\t%.2f\t%.2f\n" % (pGA, pTT, pGG, pGC))
output.write("%.2f\t%.2f\t%.2f\t%.2f\n" % (pCA, pCT, pCG, pCC))

output2 = open("Total-Allele-Counts" + args.output, 'w')
output2.write("A\tT\tG\tC\n")
output2.write("%s\t%s\t%s\t%s\n" % (totalA, totalT, totalG, totalC))

print("\nA\tT\tG\tC")
print("%s\t%s\t%s\t%s\n" % (totalA, totalT, totalG, totalC))

print("%.2f\t%.2f\t%.2f\t%.2f" % (pAA, pAT, pAG, pAC))
print("%.2f\t%.2f\t%.2f\t%.2f" % (pTA, pTT, pTG, pTC))
print("%.2f\t%.2f\t%.2f\t%.2f" % (pGA, pTT, pGG, pGC))
print("%.2f\t%.2f\t%.2f\t%.2f" % (pCA, pCT, pCG, pCC))

datafile.close()
output.close()
ref.close()
print('Done!')
