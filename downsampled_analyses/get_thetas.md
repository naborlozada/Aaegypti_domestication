
## Calculate the genetic diversity and Tajima's D

1) Calculate genotype likelihoods and the site allele frequency (SAF) 

```bash
myWDIR=/scr/core/nlozada/aedes_aegypti/scripts/;
# list with a full path for a TXT file containing the WGS BAM alignments of each inidividual in a population 
popLIST=my_aaeg_populations.list.txt;
# reference genome (in fasta format)
refGenome=/scr/core/nlozada/aedes_aegypti/reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta;

# calculate genotype likelihoods using the GATK model

# main output extension will be *.saf
for i in $popLIST; do
   angsd -b $i -anc $refGenome -out /scr/core/nlozada/aedes_aegypti/outputs/$i.angsd -ref $refGenome -minMapQ 10 -minQ 10 -minInd 1 -doSaf 1 -GL 2 -nThreads 8;
   wait;
   sleep 2;
done

wait;
```

2) Calculate the site frequency spectrum (SFS)

```bash
popSAFindex=scr/core/nlozada/aedes_aegypti/outputs/*.angsd.saf.idx;

# main output extension will be *.saf
for i in $popSAFindex; do
   realSFS  $popSAFindex  -maxiter 100  -cores 60 > $popSAFindex.saf2sfs.txt 2>&1 | tee $popSAFindex.saf2sfs.stderr.log;
   wait;
   sleep 2;
done

wait;
```
  
3) sdsdsd 
