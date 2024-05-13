## Protocol detection of Selection

All 14,677 protein coding genes of the *Aedes aegypti* (AaegL5) were screened for Selection based on the dN/dS ratio based on the YN model implemented in PAML package. This approach requires the comparison of two sister species. A summary steps to achive this goal are described below.

1) Detection of orthologues genes between *Ae. aegypti* and *Ae. albopictus* (AalboF and AalboFPA) using `proteinortho v6.3.0`. Approximately, 11,800 orthologous genes were detected between both species.
2) Protein coding sequences files were created for each single pair of orthologous genes, and were translated to aminoacid sequences using `transeq` from `EMBOSS v6.6.0.0`.
3) Codon alignments were created and refined by **removing stop codons** using `macse v2.07`. Next, each codon alignment was parsed `macse v2.07`. 
4) The YN model from PAML was run with the parameters (see below).


```bash
# 1.
# get graph matrix from reciprocal best blast hist (rBBHs):
./proteinortho -project=protein_ortho_aaegL5_vs_alboFs   -p=blastp+ -cpus=60 -sim=1 -singles -xml -identity=0.25 -coverage=50 evalue=0.00001  orthologs/aaegL5_vs_aalboX/genomes/*.faa  2>&1 | tee aedes-aegypti.protein_ortho_aaegL5_vs_alboFs.stderr.log

# Next, cluster all rBBHs using as "input" the project outfile and get orthhologs:
./proteinortho_clustering protein_ortho_aaegL5_vs_alboFs.blast-graph

# 2.
# Make file alignments
# Fasta sequences from each population VCF. SNP vcf positions were transform to protein coding sequences using the program vcf2fasta from Santiago Sanchez (https://github.com/santiagosnchez/vcf2fasta.git). 
# Fasta sequences from Ae. albopictus from two VectorBase reference genomes were used two join them with their corresponding orthologous gene.
# Merge orthologues genes using custom PERL script:



# 3.

# 4.


```

```R
.libPath();


