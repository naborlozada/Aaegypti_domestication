#!/usr/bin/perl

use strict;
use warnings;




# LOCAL SELETION
# ------------------------------------------------------------------------------------------------------------------------------------------------
# PERL script:
#      perl  genome_CDS_extraction.downsample_n_orthoGenes.v1.pl  american_samoa.tafuna_village_CDS


if ($#ARGV<0) {
    print "\n\tUsage:  perl genome_CDS_extraction.downsample_n_orthoGenes.v1.pl  POP_NAME \n";
    print "\n\n  Only one argument is mandatory: POP_NAME = Population target DIR_NAME.\n\n";
    exit (0);
}  

# ------------------------------------------------------------------------------------------------------------------------------------------------

#my $POP_TARGET = "kenya.virhembe_CDS";
my $POP_TARGET = $ARGV[0];
my ($COUNTRY,$POPULATION) = split /\./, $POP_TARGET;
   $POPULATION=~s/\_CDS//gis;

print "\n\nTARGET_POP: $POP_TARGET\n\n";





my $localdate = localtime();

my %HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_ALBOPICTUS;
my %HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_byPOP;
my %HASH_AAEGYPTI_REFGENOME_FASTASEQ; 
my %HASH_AALBOPICTUS_REFGENOME_FASTASEQ;
my %MAIN_ORTHOS_EXTRACTION;
my %HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES; 
my %hash_downsample_mosquitoes_LIST;

my %HASH_MATCH_SAMPLE_NAMES_TO_SAMPLE_ID;


my $LOC_WORKING_DIR="/home/nlozada/aedes_aegypti/NEE_paper/results/selection_dnds/populations_vcffastas_mktest/scripts/downsampling";
my $SUMMARY_OUTPUTS_DIR="/home/nlozada/aedes_aegypti/NEE_paper/results/selection_dnds/populations_vcffastas_mktest/output.summary_files";



# /// MAIN INFILES ///
my $protein_orthos_parsed_list = "/home/nlozada/aedes_aegypti/NEE_paper/results/orthologs/aedes_aegypti.orthologues_genes.ae_aegypti__vs__ae_albopictus.whole_genome_comparison.two_assemblies.table.txt";



my $all_populations_sample_ID_names = "/home/nlozada/aedes_aegypti/NEE_paper/results/populations_vcffastas_mktest/scripts/aedes_aegypti.all_populations_samples.ID_list.txt";
open ALL_POPS_SAMPLE_IDS, $all_populations_sample_ID_names or die "CANNOT open the INFILE ALL_POPS_SAMPLE_IDS ~ $!\n";


while (<ALL_POPS_SAMPLE_IDS>) {
    chomp $_;
    my $line = $_;
    if ($line!~/^#.+/ && $line!~/^POPULATION_DIR_NAME/) {
        my ($dir_name,$sample_file) = split /\t/, $line;
           $sample_file =~s/\.vcf\.gz//gis;
        my $sample_ID    = "";
        my $sample2match = "";

        if ($sample_file=~/SRR/) {  my @split_name = split /\./, $sample_file;  $sample_ID = $split_name[$#split_name];  $sample2match = join '.', @split_name[0..$#split_name-1];  }
        else {  $sample_ID = $sample_file;  $sample2match=$sample_file;  }

        $HASH_MATCH_SAMPLE_NAMES_TO_SAMPLE_ID{$sample2match} = $sample_ID;

        #print "$line\t*$sample_ID*\t--->\t$sample2match\t|$sample_file|\n";        
    }
}



# All orthologous protein coding genes between AaegL5 vs AalboF (list of gene IDs per population DIR: FINAL SET OF GENEIDS)
# All downsampled mosquitoes list
my $aaegL5_populations_samples_downsized_LIST ="/home/nlozada/aedes_aegypti/NEE_paper/files/aedes_aegypti.downsampled_populations.minsize10.info_per_sample.txt";


# Read and get downsample list of mosquitoes:
# --------------------------------------------------------------------
open DOWNSAMPLE_INDIVIDUALS_LIST, $aaegL5_populations_samples_downsized_LIST or die "CANNOT open the DOWNSAMPLE_INDIVIDUALS_LIST ~ $!\n";

while (<DOWNSAMPLE_INDIVIDUALS_LIST>) {
    chomp $_;
    my $line = $_;
    if ($line!~/^\#.+/ && $line!~/^Nr.+/ && $line!~/^\#$/) {
        my ($Nr,$Continent,$Region,$Location,$Country,$Population,$ACRN_COUNT_n_ACRN_POP,$Sample_ID) = split /\t/, $line;
		my $sample_info = "$Continent,$Region,$Location,$Country,$Population,$ACRN_COUNT_n_ACRN_POP";
		$hash_downsample_mosquitoes_LIST{$Sample_ID} = $sample_info;

        #print "|$Sample_ID|\t**$sample_info**\n";    
    }
}
close DOWNSAMPLE_INDIVIDUALS_LIST;






# Orthologous genes distribution across populations
# Read and get downsample list of mosquitoes:
# ---------------------------------------------------------------------------------
my $single_gene_byPopulation  = "aedes-aegypti.all_populations.orthologues_genes.single_geneID_by_population.list.txt";
open SINGLE_GENE_LIST_BYPOP, "$SUMMARY_OUTPUTS_DIR/$single_gene_byPopulation" or die "CANNOT open the infile SINGLE_GENE_LIST_BYPOP ~ $!\n";

while (<SINGLE_GENE_LIST_BYPOP>) {
    chomp $_;
    my $line = $_;
    if($line!~/^#.+/ && $line!~/^GENEID.+/) {
        my ($aaeg_geneID,$country,$population) = split /\t/, $line;
        #print "|$aaeg_geneID,$country,$population|\t<|$line|>\n";
        $HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_byPOP{$country}{$population}{$aaeg_geneID}++;
    }

}




# FULLSET OF ORTHOLOGOUS GENES AaegL5 and AaalbF:
# ---------------------------------------------------------------------------------
open MAIN_ORTHO_GENES_AAEG_VS_AALBO_LIST, $protein_orthos_parsed_list or die "CANNOT open the MAIN_ORTHO_GENES_AAEG_VS_AALBO_LIST ~ $!\n";

while (<MAIN_ORTHO_GENES_AAEG_VS_AALBO_LIST>) {
    chomp $_;
    my $line = $_;
    if ($line!~/^\#.+/ && $line!~/^Aedes_aegypti_ID.+/ && $line!~/^\#$/) {
        #print "***|$line|\n";
        my ($Aedes_aegypti_ID,$Aedes_albopictus_ID,$Aedes_albopictus_Assembly_associated_geneID) = split /\t/, $line;
	    my $ortho_aegypti_info = "$Aedes_albopictus_ID,$Aedes_albopictus_Assembly_associated_geneID";
	    $HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_ALBOPICTUS{$Aedes_aegypti_ID} = $ortho_aegypti_info;
    
        $MAIN_ORTHOS_EXTRACTION{$Aedes_albopictus_ID} = $Aedes_aegypti_ID;
        #print "$Aedes_aegypti_ID\t$Aedes_albopictus_ID\n";
    }
}
close MAIN_ORTHO_GENES_AAEG_VS_AALBO_LIST;






# READ AND GET CDS FROM REFERENCE GENOME AAEGL5
# ---------------------------------------------------------------------------------
my $AAEGYPTI_REF_GENOME_ALL_CDS = "/home/nlozada/aedes_aegypti/NEE_paper/programs/scripts/reference_files/reference_genome/Aedes-aegypti-LVP_AGWG_CHROMOSOMES.AaegL5_2.cds.fasta";
open AAEGYPTI_REF_GENOME_ALL_CDS_INFILE, $AAEGYPTI_REF_GENOME_ALL_CDS or die "CANNOT open the infile GENE_COUNTS_LIST_BYPOP ~ $!\n";

my $fasta_CDS_header="";

while (<AAEGYPTI_REF_GENOME_ALL_CDS_INFILE>) {
    chomp $_;
    my $line = $_;

    if ($line=~/^\>(.+)/) {
        my $tmp_line = $1;
        my $split_nr = 0; 
        my ($iso_geneID,$geneID,$gene_name,$seqID,$seqType) = ("","","","","");

        if ($tmp_line=~/name\=/) { 
            # AAEL006361-RA gene=AAEL006361 name=SCRC2 seq_id=3 type=cds
            ($iso_geneID,$geneID,$gene_name,$seqID,$seqType) = split / /, $tmp_line;
            #print "$fasta_CDS_header\n";
        }
        else {
            ($iso_geneID,$geneID,$seqID,$seqType) = split / /, $tmp_line;
            $gene_name = "name=na";
        }
        #print "$iso_geneID\t$geneID\t$gene_name\t$seqID\t$seqType\n";
        $fasta_CDS_header = $geneID;
        $fasta_CDS_header =~ s/gene\=//gis;

    }
    else {
        my $sequence = $line;
        $HASH_AAEGYPTI_REFGENOME_FASTASEQ{$fasta_CDS_header} .= $sequence;
    }

}



# CHECK CDS FROM REFERENCE GENOME AAEGL5
# ---------------------------------------------------------------------------------
# for my $aaegypti ( keys %HASH_AAEGYPTI_REFGENOME_FASTASEQ ) {
#     #print "|$aaegypti|\n"
#     #print "$aaegypti\t|$HASH_AAEGYPTI_REFGENOME_FASTASEQ{$aaegypti}|\n"
# }









# Location of orthologous genes with Ae. albopictus
# ---------------------------------------------------------------------------------------------------------------------------------------------------


my $main_DIR="";
my $main_FULL_DIR="";


my @main_directory = qw(/home/nlozada/aedes_aegypti/NEE_paper/results/orthologs/aaegL5_vs_aalboX/genomes/sequence_extraction/);

foreach my $i (@main_directory){
    chdir $i;
    #print "[1] $i ... done\n";
    opendir DIR, ".";
    
    my @subdir= readdir DIR;
    &GoDown1(@subdir);
    closedir DIR;
}


# >>> GO to NEXT SUBDIR >>>
sub GoDown1{
    my @subdir= ();
       #system "pwd";
    @subdir= @_;

    if ($#subdir > 1){
        @subdir= @subdir[0..$#subdir];
        my @subdir2= sort {lc($a) cmp lc($b)} @subdir;

        foreach my $j (@subdir2){

            if ($j!~/^\.$/ && $j!~ /^\.\.$/ && $j!~/^\.listing$/ && $j!~/readme/ && $j!~/^.+\.log$/) {
                #print ">>>$j\n";
                my $subpwd = "/home/nlozada/aedes_aegypti/NEE_paper/results/orthologs/aaegL5_vs_aalboX/genomes/sequence_extraction/$j/";
		        $main_DIR=$j;
                $main_DIR=~s/\_CDS$//gis;
                $main_FULL_DIR = $subpwd; 
                $main_FULL_DIR =~ s/\/$//g;
                chdir $subpwd;
                opendir SUBDIR, $subpwd or die "CANNOT open the SUB-DIR $subpwd - $! \n";
                
                my @files2=();
		            @files2=readdir SUBDIR;
                #print "[2] [$subpwd]\t|$j|\n";
                &GoGoDown2(@files2);
                closedir SUBDIR;
            }
        }
    }
}


# >>> GO-GO to NEXT SUBDIR >>>
sub GoGoDown2{
    my @subdir= ();
    @subdir= @_;

    if ($#subdir > 1){
        @subdir= @subdir[0..$#subdir];
        my @subdir2= sort {lc($a) cmp lc($b)} @subdir;

        foreach my $j (@subdir2){
            chomp $j;
            my ($alignment_nts,$consensus,$codon_aln,$aa_aln) = ("","");

            if ( $j!~ /^\.$/ && $j!~ /^\.\.$/ && $j!~ /^\.listing$/ && $j=~/(.+)\.longest_isoforms\.fastas\.cds\.parsed\.fasta$/ ) {
                my $index_lines=0;            
                #print "$j\t$main_FULL_DIR\n";

                open NTS_FASTA, $j or die "CANNOT open the file $j - $!\n";

                my $fasta_header;

                # MAIN
                while (<NTS_FASTA>) {
                    chomp;
                    my $line = $_;
                    $index_lines++;

                    # fasta heafer
                    if ($line=~/^>(.+)\..+/) {
                        $fasta_header = $1;
                        #print "$main_DIR\t$line\t$fasta_header\n";
                    }
                    elsif ($line=~/^>(.+)/) {
                        $fasta_header = $1;
                        $fasta_header =~s/\..+//g; 
                        #print "$main_DIR\t$line\t$fasta_header\n";
                    }
                    else {
                        my $sequence = $line;
                        $HASH_AALBOPICTUS_REFGENOME_FASTASEQ{$main_DIR}{$fasta_header} .= $sequence;
                    }
                }

            }
        }
    }
}


print "Done getting Ae. aegypti orthologous genes from Ae. albopictus...\n";








my %MATCH_aegypti_albopictus_ortho_genes;
###   $MAIN_ORTHOS_EXTRACTION{$AALBF_single_gene} = $AAEL_single_gene;

for my $strain ( keys %HASH_AALBOPICTUS_REFGENOME_FASTASEQ ) {
    for my $fastaseq ( keys %{ $HASH_AALBOPICTUS_REFGENOME_FASTASEQ{$strain} } ) {
        
        my $albopictus_gene = "$fastaseq\t$HASH_AALBOPICTUS_REFGENOME_FASTASEQ{$strain}{$fastaseq}";
        #print "$albopictus_gene\n";
        
        # find orthologs: ~12,058
        if ( $MAIN_ORTHOS_EXTRACTION{$fastaseq} ) {
            my $aegypti_gene = $MAIN_ORTHOS_EXTRACTION{$fastaseq};
            #print "$fastaseq\t$aegypti_gene\t$strain\t$albopictus_gene\n";
            
            # this HASH == HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_ALBOPICTUS
            $MATCH_aegypti_albopictus_ortho_genes{$aegypti_gene} = $albopictus_gene;
        }
    
    }
}


print "Done parsing Ae. aegypti orthologous genes from Ae. albopictus. Ready to go through each single population fasta gene sequence an extract info.\n";


#$HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_byPOP{$country}{$population}{$aaeg_geneID}++;
#$HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_ALBOPICTUS{$Aedes_aegypti_ID} = $ortho_aegypti_info;
















# read directories
# ---------------------------------------------------------------------------------------------------------------------------------------------------

my @main_directory_VCFs = qw(/home/nlozada/aedes_aegypti/NEE_paper/results/populations.vcfs.single_sample_splits/);

foreach my $i (@main_directory_VCFs){
    chdir $i;
    #print "[1] $i ... done\n";
    opendir DIR, ".";
    
    my @subdir= readdir DIR;
    &VCFGoDown1(@subdir);
    closedir DIR;
}


# >>> GO to NEXT SUBDIR >>>
sub VCFGoDown1{
    my @subdir= ();
       #system "pwd";
    @subdir= @_;

    if ($#subdir > 1){
        @subdir= @subdir[0..$#subdir];
        my @subdir2= sort {lc($a) cmp lc($b)} @subdir;

        foreach my $j (@subdir2){

            if ($j!~/^\.$/ && $j!~ /^\.\.$/ && $j!~/^\.listing$/ && $j!~/readme/ && $j!~/^.+\.log$/ && $j!~/varcalls.selection.scripts/ && $j eq $POP_TARGET ) {
		        #my ($country,$population,$chromosome);
                #$population = $j;
                my $subpwd = "/home/nlozada/aedes_aegypti/NEE_paper/results/populations.vcfs.single_sample_splits/$j/";
		            $main_DIR=$j;
                $main_FULL_DIR = $subpwd; 
                $main_FULL_DIR =~ s/\/$//g;
                chdir $subpwd;
                opendir SUBDIR, $subpwd or die "CANNOT open the SUB-DIR $subpwd - $! \n";
                
                my @files2=();
		            @files2=readdir SUBDIR;
                #print "[2] [$subpwd]\t|$j|\n";
                &VCFGoGoDown2(@files2);
                closedir SUBDIR;
            }
        }
    }
}


# >>> GO-GO to NEXT SUBDIR >>>
sub VCFGoGoDown2{
    my @subdir= ();
    @subdir= @_;

    if ($#subdir > 1){
        @subdir= @subdir[0..$#subdir];
        my @subdir2= sort {lc($a) cmp lc($b)} @subdir;

        my $POPdir = $main_DIR;

        foreach my $j (@subdir2){
            chomp $j;
            my ($codon_alignment,$target_genome_sample) = ("","");
            # target file: American23.vcf2fas.ref_genome.cds.fasta

            # each target fasta file has 14677 genes in fasta format
            # ------------------------------------------------------------------------------------------
            if ( $j!~ /^\.$/ && $j!~ /^\.\.$/ && $j!~ /^\.listing$/ && $j!~ /^\.vcf\.gz$/ && $j!~ /^\.tbi$/ && $j!~/.+\.vcf2fas\.ref_genome\.fasta$/ && $j!~/.+\.vcf2fas\.ref_genome\.fasta\.index$/ && $j=~/(.+)\.vcf2fas\.ref_genome\.cds\.fasta$/ ) {
                my $sampleID = $1;
                   $target_genome_sample = $j;
                
                
                open INFILE, "$target_genome_sample" or die "CANNOT open INFILE $target_genome_sample ~ $!\n";
                my $fasta_header = "";

                # MAIN
                # ------------------------------------------------------------------------------------------
                while (<INFILE>) {
                    chomp;
                    my $line = $_;

                    #if ($line=~/^>.+/) {
                    #    $fasta_header = "$sampleID\.$POPdir";
                    #    
                    #}
                    if ($line=~/^\>(.+)/) {
                        my $tmp_line = $1;
                        my $split_nr = 0; 
                        my ($iso_geneID,$geneID,$gene_name,$seqID,$seqType) = ("","","","","");
                        my @split_line_info = split / /, $line;
                           $geneID = $split_line_info[1];

                        $fasta_header = $geneID;
                        $fasta_header =~ s/gene\=//gis;

                        #print "$POPdir\t$sampleID\t**$j**\t$fasta_header\n";
                        
                        #if ($tmp_line=~/name\=/) { 
                        #    # AAEL006361-RA gene=AAEL006361 name=SCRC2 seq_id=3 type=cds
                        #    ($iso_geneID,$geneID,$gene_name,$seqID,$seqType) = split / /, $tmp_line;
                        #    #print "$fasta_CDS_header\n";
                        #}
                        #else {
                        #    ($iso_geneID,$geneID,$seqID,$seqType) = split / /, $tmp_line;
                        #    $gene_name = "name=na";
                        #}
                        #print "$iso_geneID\t$geneID\t$gene_name\t$seqID\t$seqType\n";
                    }
                    else {
                        my $sequence = $line;
                        $HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES{$POPdir}{$main_FULL_DIR}{$fasta_header}{$sampleID} .= $sequence;
                    }
                }
                # ------------------------------------------------------------------------------------------
            }
        }
    }
}


print "Done parsing each single gene from each population of Ae. aegypti. Now, parse consensus name.\n";






my $counter_ortho_w_snps = 0;
my $counter_ortho_no_snps = 0;
my $counter_no_ortho_snps_unknown = 0;
my %counter_CDS_single_sample; 

my $DIR_downsampled = "protein_coding_sequences.downsampled";
my $DIR_NON_downsampled = "protein_coding_sequences.non_downsampled";


my %HASH_FINAL_GROUPING_ALL_SEQUENCES;

for my $popdir ( keys %HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES ) {
    
    for my $fullDIR ( keys %{ $HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES{$popdir} } ) {

        for my $aaegL5_gene_ID ( keys %{ $HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES{$popdir}{$fullDIR} } ) {
            my $ortho_FLAG = "";


            # MAKE OUTFILE FORMAT NAME:
            # ---------------------------------------------------------------------------------------------------------
            # what genes have orthologs:
            if ( $HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_ALBOPICTUS{$aaegL5_gene_ID} ) {
                # NOW, what of these ortho-genes have SNPS:
                if ( $HASH_AAEGYPTI_GENEIDS_W_ORTHOLOGOUS_byPOP{$COUNTRY}{$POPULATION}{$aaegL5_gene_ID} ) {    
                    $counter_ortho_w_snps++;
                    $ortho_FLAG="ortho.snps";
                }
                else {
                    $counter_ortho_no_snps++;
                    $ortho_FLAG="ortho.no_snps";
                    #print ">>>|$aaegL5_gene_ID|\n";
                }
            }
            else {
                $counter_no_ortho_snps_unknown++;
                $ortho_FLAG="no_ortho.unknown";
            }
            
            my $outfile_yesDownsampled = "$aaegL5_gene_ID\.cds\.downsampled\.$ortho_FLAG\.fasta";
            my $outfile_notDownsampled = "$aaegL5_gene_ID\.cds\.non_downsampled\.$ortho_FLAG\.fasta";
            #print ">|$aaegL5_gene_ID|\t$fullDIR\t$outfile\n";
            
            # create all sequence alignments, even if they do not have an ortho ("no_ortho.unknown") but also 

            # DOWNSAMPLED:
            open OUTIFILE_FASTA_DOWNSAMPLED, ">>$fullDIR/$DIR_downsampled/$outfile_yesDownsampled" or die "CANNOT open OUTIFILE_FASTA_DOWNSAMPLED $aaegL5_gene_ID ~ $! \n";
            
            # NON-DOWNSAMPLED:
            open OUTIFILE_FASTA_NON_DOWNSAMPLED, ">>$fullDIR/$DIR_NON_downsampled/$outfile_notDownsampled" or die "CANNOT open OUTIFILE_FASTA_NON_DOWNSAMPLED $aaegL5_gene_ID ~ $! \n";
            # ---------------------------------------------------------------------------------------------------------
            #perl  genome_CDS_extraction.downsample_n_orthoGenes.v1.pl  american_samoa.tafuna_village_CDS | grep '^>'  | awk '{if($2~/cds.ortho.snps/){print "ORTHO_SNPS\t"$0} else if($2~/cds.ortho.no_snps/){print "ORTHO_NO_SNPS\t"$0} else if($2~/cds.no_ortho.unknown/){print "NO_ORTHO_UNKNOWN\t"$0} else{print "WARNING\t"$0} }' | cut -f 1 | sort | uniq -c | sort -n
            #    646 ORTHO_NO_SNPS
            #   2619 NO_ORTHO_UNKNOWN
            #  11412 ORTHO_SNPS




            # ---------------------------------------------------------------------------------------------------------
            # [STEP-1] ADD "REFERENCE GENOME" GENES
            # ---------------------------------------------------------------------------------------------------------
            if ( $HASH_AAEGYPTI_REFGENOME_FASTASEQ{$aaegL5_gene_ID} ) {
                my $aaegL5_gene_SEQ = $HASH_AAEGYPTI_REFGENOME_FASTASEQ{$aaegL5_gene_ID};
                my $refgenome_geneID = "$aaegL5_gene_ID\|REF_AAEGL5";
                #print "$popdir\t>$refgenome_geneID\n**$aaegL5_gene_SEQ**\n";
                print OUTIFILE_FASTA_DOWNSAMPLED ">$refgenome_geneID\n$aaegL5_gene_SEQ\n";
                print OUTIFILE_FASTA_NON_DOWNSAMPLED ">$refgenome_geneID\n$aaegL5_gene_SEQ\n";
            }


            for my $fasta_header_sample_ID ( keys %{ $HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES{$popdir}{$fullDIR}{$aaegL5_gene_ID} } ) {
                my $new_fasta_header_sample_ID = $fasta_header_sample_ID;
                   #$new_fasta_header_sample_ID =~ s/\./_/gis;
                   $new_fasta_header_sample_ID =~ s/jeddah\_al\-rawabi/jeddah/gis;
                   $new_fasta_header_sample_ID =~ s/\_CDS$//gis;


                my $UNIQUE_SAMPLE_ID="";
                # MAIN SHORT_UNIQ_SAMPLE_ID
                if ( $HASH_MATCH_SAMPLE_NAMES_TO_SAMPLE_ID{$new_fasta_header_sample_ID} ) {  $UNIQUE_SAMPLE_ID = $HASH_MATCH_SAMPLE_NAMES_TO_SAMPLE_ID{$new_fasta_header_sample_ID};  }
                else {  $UNIQUE_SAMPLE_ID  = $new_fasta_header_sample_ID; }
                #print ">|$aaegL5_gene_ID|\t$fasta_header_sample_ID\t*$UNIQUE_SAMPLE_ID*\n";


                #my $codon_sequences = ">$new_fasta_header_sample_ID\|$HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES{$popdir}{$fullDIR}{$aaegL5_gene_ID}{$fasta_header_sample_ID}";

                #print "$popdir\t$aaegL5_gene_ID\t$fasta_header_sample_ID\t**$codon_sequences**\n";
                #$HASH_FINAL_GROUPING_ALL_SEQUENCES{$POP_TARGET}{$aaegL5_gene_ID} = $codon_sequences;
                $counter_CDS_single_sample{$UNIQUE_SAMPLE_ID}++;


                # ---------------------------------------------------------------------------------------------------------
                # [STEP-2] ADD "SINGLE MOSQUITOES" GENES
                # ---------------------------------------------------------------------------------------------------------
                # ALL SAMPLES (non-downsamples)
                my $codon_sequences = $HASH_AAEGYPTI_SINGLE_MOSQUITO_GENES{$popdir}{$fullDIR}{$aaegL5_gene_ID}{$fasta_header_sample_ID};
                my $fasta_header = "$aaegL5_gene_ID\|$UNIQUE_SAMPLE_ID";
                print OUTIFILE_FASTA_NON_DOWNSAMPLED ">$fasta_header\n$codon_sequences\n";
                #print  ">$UNIQUE_SAMPLE_ID\t$fasta_header\t**$codon_sequences**\n";


                # ALL SAMPLES (downsamples)
                if ( $hash_downsample_mosquitoes_LIST{$UNIQUE_SAMPLE_ID} ) {
                    #print ">>>$popdir\t$UNIQUE_SAMPLE_ID\t$aaegL5_gene_ID\t$fasta_header\t**$codon_sequences**\n";
                    print OUTIFILE_FASTA_DOWNSAMPLED ">$fasta_header\n$codon_sequences\n";
                }
                elsif ( $POP_TARGET eq "american_samoa.tafuna_village_CDS" && $UNIQUE_SAMPLE_ID=~/American16/) {
                    print OUTIFILE_FASTA_DOWNSAMPLED ">$fasta_header\n$codon_sequences\n";
                }
                elsif ( $POP_TARGET eq "senegal.ngoye_CDS" && $UNIQUE_SAMPLE_ID=~/SRR11006762/) {
                    # ref.downsampled.ID SRR11006760 not present (possibly due this sample has NO SNPs in their sequences or so?), so I picked another one: XXXXX
                    print OUTIFILE_FASTA_DOWNSAMPLED ">$fasta_header\n$codon_sequences\n";
                }
                else {
                    #print ">>>$popdir\t$aaegL5_gene_ID\t*$UNIQUE_SAMPLE_ID*\t$fasta_header_sample_ID\n";
                }

            }


            # ---------------------------------------------------------------------------------------------------------
            # [STEP-3] ADD "REFERENCE GENOME" GENES
            # ---------------------------------------------------------------------------------------------------------
            if ( $MATCH_aegypti_albopictus_ortho_genes{$aaegL5_gene_ID} ) {
                my ($aalboF_gene_ID,$aalboF_gene_SEQ) = split /\t/, $MATCH_aegypti_albopictus_ortho_genes{$aaegL5_gene_ID};
                my $fasta_header = "$aaegL5_gene_ID\|$aalboF_gene_ID";
                #print "$popdir\t>$refgenome_geneID\n**$aaegL5_gene_SEQ**\n";
                print OUTIFILE_FASTA_DOWNSAMPLED ">$fasta_header\n$aalboF_gene_SEQ\n";
                print OUTIFILE_FASTA_NON_DOWNSAMPLED ">$fasta_header\n$aalboF_gene_SEQ\n";
            }

            close OUTIFILE_FASTA_DOWNSAMPLED;
            close OUTIFILE_FASTA_NON_DOWNSAMPLED;
        }
    }
}


print "Done parsing populations information of Ae. aegypti. Now, make output sequence fasta files..\n";
print "\n\n";


