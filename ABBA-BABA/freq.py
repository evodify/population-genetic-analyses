#!/usr/bin/env python2
"""
Original script was written by Simon Martin (shm45@cam.ac.uk)
Publication: MARTIN ET AL. 2013 GENOME-WIDE EVIDENCE FOR SPECIATION WITH GENE FLOW IN HELICONIUS BUTTERFLIES
Email: shm45@cam.ac.uk
07-04-2013

This script calculates the frequency of either the minor allele or the derived allele in each defined population.
To calculate the derived allele frequency, an out-group sequence/population must be defined as a reference for the ancestral state.
If multiple out-groups exist and they are not fixed for a single allele, a consensus can be taken.
Only biallelic sites will be considered.

Input "calls file" format:
scaffold  position  ind1  ind2  ind3  etc...
scf1  1 T T C Y

NOTE a header row with unique names for all individuals is essential.

Output format is .csv with scaffold, position and frequencies. "NA" will be specified for invalid or monomorphic sites.

Command:

python freq.py -i <input file> -o <output file> -p <population_string> -a derived -O <out-group population name> --consensus -M 3

the populations are specified in a single string as follows:

-p "pop1Name[ind1,ind2,ind3,ind4];pop2Name[ind1,ind2,ind3,ind4]" (quotation marks must be present to avoid conflict with unix)

the pout-group must match one of the population names. e.g. -O pop2Name

***********************************************************************************************

Modified to use single nucleotide gaps "-" by Dmytro Kryvokhyzha (dmytro.kryvokhyzha@evobio.eu)
"""

import sys

### Functions
def get_intv(string,borders = "()",inc = False):
  if len(borders) != 2:
    print "WARNING: borders must contain two characters"
  starts = []
  ends = []
  output = []
  for x in range(len(string)):
    if string[x] == borders[0]:
      starts.append(x)
    if string[x] == borders[1]:
      ends.append(x+1)
  if len(starts) <= len(ends):
    for n in range(len(starts)):
      if inc:
        output.append(string[starts[n]:ends[n]])
      else:
        output.append(string[starts[n]+1:ends[n]-1])
  else:
    for n in range(len(ends)):
      if inc:
        output.append(string[starts[n]:ends[n]])
      else:
        output.append(string[starts[n]+1:ends[n]-1])
  return output

def haplo(calls):
  output = []
  for call in calls:
    if call in "ACGTN-":
      output.append(call)
      output.append(call)
    elif call == "K":
      output.append("G")
      output.append("T")
    elif call == "M":
      output.append("A")
      output.append("C")
    elif call == "R":
      output.append("A")
      output.append("G")
    elif call == "S":
      output.append("C")
      output.append("G")
    elif call == "W":
      output.append("A")
      output.append("T")
    elif call == "Y":
      output.append("C")
      output.append("T")
    else:
      print "WARNING", call, "is not recognised as a valid base or ambiguous base"
      output.append("N")
      output.append("N")
  return output

def getOptionValue(option):
  optionPos = [i for i,j in enumerate(sys.argv) if j == option][0]
  optionValue = sys.argv[optionPos + 1]
  return optionValue
  
def unique(things):
  output = []
  for x in things:
    if x not in output:
      output.append(x)
  output.sort()
  return output


def exclude(things, x):
  output = [i for i in things if i != x]
  return(output)

def uniqueAlleles(bases):
  haploBases = haplo(bases)
  output = unique([i for i in haploBases if i in "ACGT-"])
  return output


def baseFreq(bases,base):
  haploBases = [i for i in haplo(bases) if i in "ACGT-"]
  freq = (float(haploBases.count(base))) / len(haploBases)
  return freq

def mostCommon(things):
  output = []
  counts = []
  uniqueThings = unique(things)
  for thing in uniqueThings:
    counts.append(things.count(thing))
  maxCount = max(counts)
  for n in range(len(counts)):
    if counts[n] == maxCount:
      output.append(uniqueThings[n])
  return output

def majorAllele(bases):
  haploBases = [i for i in haplo(bases) if i in "ACGT-"]
  major = mostCommon(haploBases)[0]
  return major[0]

### get files

if "-i" in sys.argv:
  fileName = getOptionValue("-i")
else:
  print "\nplease specify input file name using -i <file_name> \n"
  sys.exit()
  
file = open(fileName, "rU")
line = file.readline()
names = line.split()
line= file.readline()


if "-o" in sys.argv:
  outName = getOptionValue("-o")
  out = open(outName, "w")
else:
  print "\nplease specify output file name using -o <file_name> \n"
  sys.exit()

if "-p" in sys.argv:
  popString = getOptionValue("-p")
else:
  print "\nplease specify populations using -p\n"
  sys.exit()

if "-a" in sys.argv:
  if getOptionValue("-a") == "derived":
    derived = True
  elif getOptionValue("-a") == "minor":
    derived = False
  else:
    print "\nAllele (-a) can only be 'minor' or 'derived'."
    sys.exit()
else:
  print "\nPlease specify whether the derived or minor allele frequency must be calculated.\n"
  sys.exit()

if "-O" in sys.argv:
  outGroup = getOptionValue("-O")
else:
  if derived:
    print "\nPlease specify outgroup population using -O\n"
    sys.exit()

if "-M" in sys.argv:
  popMin = int(getOptionValue("-M"))
else:
  popMin = 1


if "--consensus" in sys.argv:
  outgroupConsensus = True
else:
  outgroupConsensus = False

pops = []
#for each population, store the name and individual names
for popData in popString.strip("\"").split(";"):
  currentPop = popData.split("[")[0]
  pops.append(currentPop)
  vars()[currentPop + "Inds"] = get_intv(popData,"[]")[0].split(",")
  for ind in vars()[currentPop + "Inds"]:
    if ind not in names:
      print ind, "not found in header line."
      sys.exit()

if derived and outGroup not in pops:
  print "\nThe specified outgroup, ", outGroup, ", was not a specified population."
  sys.exit()

# write output header
out.write(names[0] + "," + names[1])
for pop in pops:
  out.write("," + pop)
out.write("\n")

linesDone = 0

### for each line, check if its a biallelic SNP, if so, continue to other populations
while len(line) > 1:
  objects = line.split()
  output = [objects[0],objects[1]]
  # check not triallelic
  alleles = uniqueAlleles(objects[2:])
  if len(alleles) == 0 or len(alleles) > 2:
    for pop in pops:
      output.append("NA")
  elif len(alleles) == 1:
    for pop in pops:
      output.append("0.0")
  else:
    # get major allele or ancestral state
    if derived:
      ogBases = []
      for ind in vars()[outGroup + "Inds"]:
        ogBases.append(objects[names.index(ind)])
      ogAlleles = uniqueAlleles(ogBases)
      if len(ogAlleles) == 1:
        refState = ogAlleles[0]
      elif len(ogAlleles) == 2 and outgroupConsensus:
        refState = majorAllele(ogBases)
      else:
        refState = None
    else:
      refState = majorAllele(objects[2:])
    if refState:
      for pop in pops:
        popCalls = []
        for ind in vars()[pop + "Inds"]:
          popCalls.append(objects[names.index(ind)])
        if len(exclude(popCalls,"N")) >= popMin:
          freq = 1 - baseFreq(popCalls,refState)
        else: freq = "NA"
        output.append(str(freq))
    else:
      for pop in pops:
        output.append("NA")
  out.write(",".join(output))
  out.write("\n")    
  line = file.readline()
  linesDone += 1
  if linesDone % 1000000 == 0:
    print linesDone, "lines done..."

out.close
file.close
