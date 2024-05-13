## Protocol detection of Selection

All 14,677 protein coding genes of the *Aedes aegypti* (AaegL5) were screened for Selection based on the dN/dS ratio based on the YN model implemented in PAML package. This approach requires the comparison of two sister species. A summary steps to achive this goal are described below.

1) Detection of orthologues genes between *Ae. aegypti* and *Ae. albopictus* (AalboF and AalboFPA) using `proteinortho v6.3.0`. Approximately, 11,800 orthologous genes were detected between both species.
2) Protein coding sequences files were created for each single pair of orthologous genes, and were translated to aminoacid sequences using `transeq` from `EMBOSS v6.6.0.0`.
3) Codon alignments were created and refined by **removing stop codons** using `macse v2.07`. Next, each codon alignment was parsed `macse v2.07`. Final refinement of the alignment was performed using "pal2nal.pl"
4) The YN model from PAML was run with the parameters (see below).


```bash
# 1. Get orthologues genes:
# get graph matrix from reciprocal best blast hist (rBBHs):
./proteinortho -project=protein_ortho_aaegL5_vs_alboFs   -p=blastp+ -cpus=60 -sim=1 -singles -xml -identity=0.25 -coverage=50 evalue=0.00001  orthologs/aaegL5_vs_aalboX/genomes/*.faa  2>&1 | tee aedes-aegypti.protein_ortho_aaegL5_vs_alboFs.stderr.log

# Next, cluster all rBBHs using as "input" the project outfile and get orthhologs:
./proteinortho_clustering protein_ortho_aaegL5_vs_alboFs.blast-graph

# 2. Make file alignments:
# Fasta sequences from each population VCF. SNP vcf positions were transform to protein coding sequences using the program vcf2fasta from Santiago Sanchez (https://github.com/santiagosnchez/vcf2fasta.git). 
# Fasta sequences from Ae. albopictus from two VectorBase reference genomes were used two join them with their corresponding orthologous gene.
# Merge orthologues genes using custom PERL script:
perl merge_orthologues_genes.pl

# 3. Codon alignment, refinament and parse final output:
# make alignment
java -jar macse_v2.07.jar -prog  alignSequences  -seq  POPNAME.AAEL000001.orthos.fas  -out_NT POPNAME.AAEL000001.orthos_NT.fas  -out_AA POPNAME.AAEL000001.orthos_AA.fas -gc_def 1 2>> POPNAME.AAEL000001.orthos_NT.stderr.log

# refine alignment
java -jar macse_v2.07.jar -prog  refineAlignment  -align  POPNAME.AAEL000001.orthos_NT.fas  -optim 1  -max_refine_iter 1  -fs 40.0  -gap_ext 3.0  -stop 60.0  -alphabet_AA Dayhoff_6  -gc_def 1  -out_AA POPNAME.AAEL000001.orthos_AA_aln.fas  -out_NT POPNAME.AAEL000001.orthos_NT.fas  2>> POPNAME.AAEL000001.orthos_NT_NT2realn.stderr.log &> POPNAME.AAEL000001.orthos_NT_NT2realn.stderr.log

# remove/mask stop codons:
java -jar macse_v2.07.jar -prog  exportAlignment  -align  POPNAME.AAEL000001.orthos_NT_aln.fas  -codonForInternalStop NNN  -codonForInternalFS ---  -charForRemainingFS ---  -out_AA POPNAME.AAEL000001.orthos_AA_aln_noFS.fas  -out_NT POPNAME.AAEL000001.orthos_NT_aln_noFS.fas  2>> POPNAME.AAEL000001.stderr.log &> POPNAME.AAEL000001.orthos_NT_aln_noFS.stderr.log

# Final parsing of the alignment:
perl pal2nal.pl  POPNAME.AAEL000001.orthos_AA_aln_noFS.fas  POPNAME.AAEL000001.orthos_NT_aln_noFS.fas  -nomismatch -nogap -output fasta 2>> POPNAME.AAEL000001.orthos_NT_aln_noFS2codon.stderr.log > POPNAME.AAEL000001.orthos_NT_aln_noFS.codon.fas

# 4. Perform the pN/pS ratio
# Using base control file from PAML GENEID_yn00.ctl in this directory:
https://github.com/naborlozada/Aaegypti_domestication/blob/main/selection/GENEID_yn00.ctl

# This file was created separately for each protein coding in each population (>500,000 control files).
```
