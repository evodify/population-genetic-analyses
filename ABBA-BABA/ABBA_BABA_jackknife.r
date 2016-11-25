#!/usr/bin/R

### Script to perform the ABBA BABA tests with a block-jacknife procedure to estimate variance
### used for MARTIN ET AL. 2013 GENOME WIDE EVIDENCE FOR SPECIATION WITH GENE FLOW IN HELICONIUS BUTTERFLIES
### Modified for the usage on Capsella data by Dmytro Kryvokhyzha (dmytro.kryvokhyzha@evobio.eu)

# Input is a csv of allele frequencies generated with freq.py script


ff <- commandArgs(trailingOnly = T)

filename <- ff[1]
testfile <- ff[2]

tests <- read.csv(testfile)
freq_table <- read.csv(filename)

attach(freq_table)

# unique(X.CHROM)
chromNames <-  as.character(unique(X.CHROM))

#create output variables
D <- numeric(length = length(tests[,1])) # D-statistic
D_err <- numeric(length = length(tests[,1])) # D std error
D_Z <- numeric(length = length(tests[,1])) # D Z-score
D_p <- numeric(length = length(tests[,1])) # D p-value
f <- numeric(length = length(tests[,1])) # f-statistic
f_err <- numeric(length = length(tests[,1])) # D std error

### normal test

#create a list of positions that will be considered during each jackknife replicate
indices <- list()
n <- 1
for (chrom in chromNames) {
  chrom_end <- max(freq_table[which(freq_table$X.CHROM == chrom),"POS"])
  starts <- seq(1,chrom_end,1000000)
  ends <- starts + 999999
  for (x in 1:length(starts)) {
    indices[[n]] <- which(freq_table$X.CHROM %in% chromNames & !(freq_table$X.CHROM == chrom & freq_table$POS >= starts[x] & freq_table$POS <= ends[x]))
    print(n)
    n <- n + 1
    }
  }

#calculate D and f statistics
for (z in 1:length(tests[,1])) {
  donor_name <- as.character(tests$donor[z])
  recipient_name <- as.character(tests$recipient[z])
  non_recipient_name <- as.character(tests$non_recipient[z])

  ABBA <- ((1 - get(non_recipient_name)) * get(recipient_name) * get(donor_name))
  BABA <- (get(non_recipient_name) * (1 - get(recipient_name)) * get(donor_name))

  maxABBA <- ((1 - get(non_recipient_name)) * get(donor_name) * get(donor_name))
  maxBABA <- (get(non_recipient_name) * (1 - get(donor_name)) * get(donor_name))

  D[z] <- (sum(ABBA, na.rm = T) - sum(BABA, na.rm = T))/(sum(ABBA, na.rm = T) + sum(BABA, na.rm = T))
  f[z] <- (sum(ABBA, na.rm = T) - sum(BABA, na.rm = T))/(sum(maxABBA, na.rm = T) - sum(maxBABA, na.rm = T))

  # jacknife to get errors for D and f
  D_pseudo <- vector(length = length(indices))
  f_pseudo <- vector(length = length(indices))
  n <- 1
  for (x in 1:length(indices)) {
    jacked_positions <- indices[[x]]
    jacked_D <- (sum(ABBA[jacked_positions], na.rm = T) - sum(BABA[jacked_positions], na.rm = T))/(sum(ABBA[jacked_positions], na.rm = T) + sum(BABA[jacked_positions], na.rm = T))
    jacked_f <- (sum(ABBA[jacked_positions], na.rm = T) - sum(BABA[jacked_positions], na.rm = T))/(sum(maxABBA[jacked_positions], na.rm = T) - sum(maxBABA[jacked_positions], na.rm = T))

    D_pseudo[n] <- D[z]*length(indices) - jacked_D*(length(indices)-1)
    f_pseudo[n] <- f[z]*length(indices) - jacked_f*(length(indices)-1)
    
    print(n)
    n <- n + 1
    }
  
  D_err[z] <- sqrt(var(D_pseudo)/length(D_pseudo))
  D_Z[z] <- D[z] / D_err[z]
  D_p[z] <- 2*pnorm(-abs(D_Z[z]))
  f_err[z] <- sqrt(var(f_pseudo)/length(f_pseudo))
  

  #then put together the results and write
  output <- cbind(tests,D,D_err,D_Z,D_p,f,f_err)

  write.csv(output, file = paste(filename, "ABBA_BABA.csv", sep="_"), row.names = F, quote = F)
  }
