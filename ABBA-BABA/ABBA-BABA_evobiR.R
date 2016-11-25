#!/usr/bin/R

library('evobiR')

ff <- commandArgs(trailingOnly = T)

for (file in ff) {
  CalcD(alignment = file, sig.test = "N")
  # CalcD(alignment = file, sig.test = "J", block.size = 100000, replicate = 100)
  message(paste('\t', file, sep=''))
  }

# CalcD('filename.fasta', sig.test = "J", block.size = 100000, replicate = 100)
