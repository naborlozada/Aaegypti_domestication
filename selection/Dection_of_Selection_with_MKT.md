## Detection of selection using the "McDonald and Kreitman test (MKT)"

1. We first performed a search of orthologs genes between Aedes aegypti and Aedes albopictus across the whole set of 14,677 genes in the reference genome AaegL5 as described in Methods.
2. Next, for each single individual sample in an African and Out of Africa population, we created a whole genome sequences using the reference genome AaegL5, and reference nucleotides were replaced by the corresponding SNPs nucleotide change using a .
3. After downsampling each population (as described in Methods), we extracted all orthologs genes in each single individual genome and population.
4. Codon alignments were created and refined using 'macse' v2.07, and verified with 'pal2nal' program.
5. We used the 'iMKT' R package and scripts provided by the [Murga-Moreno etal 2019: iMTK group paper](https://academic.oup.com/nar/article/47/W1/W283/5488529?login=false) to detect signals of Selection on each single gene alignment in each population of Africa and Out of Africa.

```bash
## /// Step 1 ///

# Blast + clustering similarities:
proteinortho_grab_proteins.pl -tofiles protein_ortho_aaegL5_vs_alboFs.proteinortho.tsv  'Aedes-aegypti-LVP_AGWG_AaegL5_2.longest_isoforms.faa'  'VectorBase-55_AalbopictusFoshanFPA.longest_isoforms.faa'  'VectorBase-61_AalbopictusFoshan.longest_isoforms.faa'  -p=blastp+  -cpus=60  -sim=1  -18 singles  -xml  -identity=0.25  -coverage=50  -evalue=0.00001

# Orthology
proteinortho_clustering  protein_ortho_aaegL5_vs_alboFs.blast-graph  >  protein_ortho_aaegL5_vs_alboFs.proteinortho-graph.main_output.txt


## /// Step 2 & 3 ///
```bash
$CDS_POPS=/home/nlozada/aaegypti/all_downsampled_pops/alignments/*_CDS

for i in CDS_POP;
   pop_fasta_genome = ${i}.fasta;
   pop_VCF_name = ${i}.vcf.gz;
   bcftools consensus  --fasta-ref Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.fasta  --output $pop_fasta_genome  $pop_VCF_name;
   agat_sp_extract_sequences.pl -g Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.GFF3 -f $pop_fasta_genome -type cds
done


## /// Step 4 ///
my $

for i in  CDS_POP;
   cd $i;
   infile.cds.fasta = $i;
   for j in infile.cds.fasta;
      java -jar macse.jar $j  



## /// Step 5 ///
