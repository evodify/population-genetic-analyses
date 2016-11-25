# ABBA-BABA test for introgression


The description of the test is provided in [Durand et al. 2011](http://mbe.oxfordjournals.org/content/28/8/2239.short).

The analysis can be performed in two ways:

## On a set of samples

The scripts are taken from [Martin et al. (2013)](http://datadryad.org/resource/doi:10.5061/dryad.dk712) with modifications.

### Generate derived allele frequency file

```
python freq.py -i genotypeCallsFile.tab -o genotypeCallsFile.freq -p "groupA[sample1,sample2,sample3];groupB[sample4,sample5,sample6];groupC[sample7,sample8,sample9];groupD[sample10];groupG[sample11,sample12]" -a derived -O groupD
```

See [freq.py](population-genetic-analyses/ABBA-BABA/freq.py) for more details.


### Run ABBA-BABA

```
Rscript ABBA_BABA_jackknife.r genotypeCallsFile.freq tests.csv 
```

where `genotypeCallsFile.freq` is the file obtained in the previous step. `test.csv` is a file defining tests to perform.
  
An example of `test.csv`:
```
non_recipient   recipient   donor
groupA	groupB	groupG
groupA	groupB	groupC

```

## On separate sequences


[ABBA-BABA_evobiR.R](ABBA-BABA_evobiR.R) performs analysis on single sequence basis in a fasta format using [evobiR](https://cran.r-project.org/web/packages/evobiR/index.html).

The script can be run on many fasta files:

```
Rscript ABBA-BABA_evobiR.R file1.fasta file2.fasta file3.fasta
```
[example.fasta](example.fasta) - an example fasta file. Note, missing data is not allowed.

**Block-jackknifing** can be adjusted in the script with options `block.size` and `replicate`. For more options, run `?CalcD` within `R`.