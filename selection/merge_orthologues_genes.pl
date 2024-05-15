#!/usr/bin/perl


# ------------------------------------------------------------------------------------------------------------------------------------------------
# Description:
#     PERL script that reads one Population directory (among many) and foreach fasta file havind a protein coding gene from Ae. aegypti, finds its ortholog and 
# put them together in a single fasta file. This process is repeated across all genes from Ae. aegypti havind orthologues gene (~12,000) detected in a TARGET POPULATION.
#
#    Note that the protein coding genes from Ae. aegypti were first extracted from SNPs in VCF files of each population using the python program "vcf2fasta" from Santiago Sanches
# (https://github.com/santiagosnchez/vcf2fasta.git) (SEQUENCES_DB1). Then, the protein coding genes from Ae. albopictus were extracted from the reference genomes from VectorBase
# (aseemblies AalbF and AalbFPA version 51 and 61, respectively)  (SEQUENCES_DB2). Overall, this PERL script find and merge orthologs searching between datasets SEQUENCES_DB1 and
#  SEQUENCES_DB2.
#
# The program requires only argument: Population target name.
# ------------------------------------------------------------------------------------------------------------------------------------------------

use strict;
use warnings;






# ARG INPUT FILE
# ------------------------------------------------------------------------------------------------------------------------------------------------
# MY TARGET
if ($#ARGV<0) {
    print "\n\tUsage:  merge_orthologues_genes.pl  POP_NAME \n";
    print "\nOnly one argument is mandatory: POP_NAME = Population target name. Example:\n\tperl merge_orthologues_genes.pl  [POP_NAME]\n\n";
    exit (0);
}

# example
#my $POP_TARGET = "kenya.virhembe_CDS";
my $POP_TARGET = $ARGV[0];
print "\n\nTARGET_POP: $POP_TARGET\n\n";

# ------------------------------------------------------------------------------------------------------------------------------------------------










# Variables
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my $localdate      = localtime();
my $population_DIR = "";    
my $main_DIR       = "";
my $main_FULL_DIR  = "";

my $main_OUTPUT_DIR = "/path/to/aedes_aegypti/results/selection_dnds/aedes_aegypti_local_adaptation.yn_count_model";

my %MAJOR_POPULATION_GROUPS;
my %aegypti_consensus_sequences;






# orthologs
# ===================================================================================================================================================

my %orthologs_aalbF_aalbFPA;
my %MAIN_ORTHOS_AEGYPTI_ALBOPICTUS;
my %MAIN_ORTHOS_EXTRACTION;

my %albopictus_strain_FASTASEQ;


# get list of orthologous genes
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my $orthologous_genes = "/path/to/proteinortho/aedes_aegypti/output/selection.GLOBAL.pairwise/protein_ortho_aaegL5_vs_alboFs.proteinortho-graph.MAIN.output.txt";


open my $INFILE, '<', $orthologous_genes;
my @indexes;

while (my $line = <$INFILE>) {
    chomp $line;

    if ($line!~/^\#.+/) {
        my @split_line = split /\t/, $line;
        my $aaegL5  = $split_line[3];
        my $aalbFPA = $split_line[4];
        my $aalbF   = $split_line[5];

        my $target_ortho="";

        # orthologs
        if ($aaegL5=~/^AAEL.+/) {
            if ($aalbFPA=~/^AALFPA\_\d+/) {
                if ($aalbFPA=~/\,/) {   my @split = split /\,/, $aalbFPA;  $target_ortho = $split[0]; } else {  $target_ortho = $aalbFPA;   }

                # push $AALFPA->@AAEL 
                $orthologs_aalbF_aalbFPA{$aaegL5} = $target_ortho;
                #print "**$target_ortho**\t$line\n";
            }
            elsif ($aalbF=~/^AALF\d+/) {
                # push $AALF->@AAEL
                if ($aalbF=~/\,/) {   my @split = split /\,/, $aalbF;  $target_ortho = $split[0]; } else {  $target_ortho = $aalbF;   }
                #print "**$target_ortho**\t$line\n";
                $orthologs_aalbF_aalbFPA{$aaegL5} = $target_ortho;
            }
        }
    }
}



# parse aegypti genes to orthos

for my $AAEL_gene (keys %orthologs_aalbF_aalbFPA) {

    my $AAEL_single_gene = "";
    my $AALBF_single_gene = "";

    if ($AAEL_gene!~/\,/) {
        $AAEL_single_gene = $AAEL_gene;
        $AALBF_single_gene = $orthologs_aalbF_aalbFPA{$AAEL_gene};
        $MAIN_ORTHOS_AEGYPTI_ALBOPICTUS{$AAEL_single_gene} = $AALBF_single_gene;
        $MAIN_ORTHOS_EXTRACTION{$AALBF_single_gene} = $AAEL_single_gene;
    } 
    else {  
        my @split = split /\,/, $AAEL_gene;
        
        foreach my $gene (@split) {
            chomp $gene;
            $AAEL_single_gene = $gene;
            $AALBF_single_gene = $orthologs_aalbF_aalbFPA{$AAEL_gene};
            #print "$AAEL_single_gene\t$orthologs_aalbF_aalbFPA{$AAEL_gene}\t*$AAEL_gene*\n";
            $MAIN_ORTHOS_AEGYPTI_ALBOPICTUS{$AAEL_single_gene} = $AALBF_single_gene;
            $MAIN_ORTHOS_EXTRACTION{$AALBF_single_gene} = $AAEL_single_gene;
        }
    }
}






# get orthologous genes
# ---------------------------------------------------------------------------------------------------------------------------------------------------

my @main_directory = qw(/path/to/output/orthologs/aaegL5_vs_aalboX/genomes/sequence_extraction/);

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
                my $subpwd = "/path/to/output/orthologs/aaegL5_vs_aalboX/genomes/sequence_extraction/$j/";
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
                        #print "$line\t$fasta_header\n";
                    }
                    elsif ($line=~/^>(.+)/) {
                        $fasta_header = $1;
                        $fasta_header =~s/\..+//g; 
                        #print "$line\t$fasta_header\n";
                    }
                    else {
                        my $sequence = $line;
                        $albopictus_strain_FASTASEQ{$main_DIR}{$fasta_header} .= $sequence;
                    }
                }

            }
        }
    }
}


print "Done getting Ae. aegypti orthologous genes from Ae. albopictus...\n";




my %MATCH_aegypti_albopictus_ortho_genes;

for my $strain ( keys %albopictus_strain_FASTASEQ ) {
    for my $fastaseq ( keys %{ $albopictus_strain_FASTASEQ{$strain} } ) {
        
        my $albopictus_gene = "$strain\|>$fastaseq\t$albopictus_strain_FASTASEQ{$strain}{$fastaseq}";
        
        # find orthologs: ~12,058
        if ( $MAIN_ORTHOS_EXTRACTION{$fastaseq} ) {
            my $aegypti_gene = $MAIN_ORTHOS_EXTRACTION{$fastaseq};
            $MATCH_aegypti_albopictus_ortho_genes{$aegypti_gene} = $albopictus_gene;
        }
    }
}


print "Done parsing Ae. aegypti orthologous genes from Ae. albopictus. Ready to go through each single population fasta gene sequence an extract info.\n";















# read main directory with all data separately by populations. 
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my @main_directory_VCFs = qw(/path/to/results/selection_dnds/aedes_aegypti_local_adaptation.yn_count_model/);

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
    @subdir= @_;

    if ($#subdir > 1){
        @subdir= @subdir[0..$#subdir];
        my @subdir2= sort {lc($a) cmp lc($b)} @subdir;

        foreach my $j (@subdir2){

            if ($j!~/^\.$/ && $j!~ /^\.\.$/ && $j!~/^\.listing$/ && $j!~/readme/ && $j!~/^.+\.log$/ && $j eq $POP_TARGET ) {
		my ($country,$population,$chromosome);
                #$country    = $1;
                $population = $j;
                #$chromosome = $3;
                my $subpwd = "/path/to/results/selection_dnds/aedes_aegypti_local_adaptation.yn_count_model/$j/";
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
        #$POPdir =~ s/\_CDS$//g;
	
	
        foreach my $j (@subdir2){
            chomp $j;
            my ($codon_alignment,$consensus,$codon_aln,$aa_aln) = ("","","","");

            # single consensus for a protein coding gene of a single population (~12,000 genes)
            # ------------------------------------
            if ( $j!~ /^\.$/ && $j!~ /^\.\.$/ && $j!~ /^\.listing$/ && $j!~/.+\.fasta$/ && $j!~/(AAEL\d+)\.fas$/ && $j=~/(.+)\.cons\.fas$/ ) {
                $consensus = $j;
                
                my $geneID = $consensus;
                   $geneID =~s/\.cons\.fas$//g;

                open INFILE, "$consensus" or die "CANNOT open INFILE $consensus ~ $!\n";
                my $fasta_header = "";

                # MAIN
                # ------------------------------------------------------------------------------------------
                while (<INFILE>) {
                    chomp;
                    my $line = $_;

                    if ($line=~/^>.+/) {
                        $fasta_header = "$geneID\.$POPdir";
                    }
                    else {
                        my $sequence = $line;
                        $aegypti_consensus_sequences{$POPdir}{$main_FULL_DIR}{$geneID}{$fasta_header} .= $sequence;
                    }
                }
            }
        }
    }
}


print "Done parsing each single gene from each population of Ae. aegypti. Now, parse consensus name.\n";









my %HASH_FINAL_GROUPING_ALL_SEQUENCES;

for my $popdir ( keys %aegypti_consensus_sequences ) {
    for my $fullDIR ( keys %{ $aegypti_consensus_sequences{$popdir} } ) {
        for my $geneID ( keys %{ $aegypti_consensus_sequences{$popdir}{$fullDIR} } ) {
            for my $fasta_header ( keys %{ $aegypti_consensus_sequences{$popdir}{$fullDIR}{$geneID} } ) {
                my $new_fasta_header = $fasta_header;
                   $new_fasta_header =~ s/\./_/gis;
                   $new_fasta_header =~ s/jeddah\_al\-rawabi/jeddah/gis;
                   $new_fasta_header =~ s/\_CDS$//gis;
                
 
                my $codon_sequences = ">$new_fasta_header\|$aegypti_consensus_sequences{$popdir}{$fullDIR}{$geneID}{$fasta_header}";

                #print "FASTA_SEQ $popdir\t$geneID\t$codon_sequences\n";
                $HASH_FINAL_GROUPING_ALL_SEQUENCES{$POP_TARGET}{$geneID} = $codon_sequences;
            }
        }
    }
}


print "Done parsing populations information of Ae. aegypti. Now, make output sequence fasta files..\n";



for my $GROUP (keys %HASH_FINAL_GROUPING_ALL_SEQUENCES) {
    my $filename = $GROUP;
       $filename =~ s/\_CDS$//gis;

    for my $GENEID (keys %{ $HASH_FINAL_GROUPING_ALL_SEQUENCES{$GROUP} } ) {

        my $outfile = ""; #"aedes_aegypti\.local\.$GROUP\.$GENEID\.orthos\.fasta";
        my $main_popdir = "";

        ## >>> FILTER!! <<<
        ## --------------------------------------------------------------------------------
        if ( $MATCH_aegypti_albopictus_ortho_genes{$GENEID} ) {    
            # grouping
            $outfile = "$filename\.$GENEID\.orthos\.fas";
            #system "rm $outfile";
            open TARGET_POPULATION, ">>$main_OUTPUT_DIR\/$POP_TARGET\/$outfile" or die "CANNOT open the outfile FINAL OUTFILE at $POP_TARGET DIR ~ $!\n";
        }
        
        # parse fasta sequences
        my $fasta_sequence_AAEGYPTI = $HASH_FINAL_GROUPING_ALL_SEQUENCES{$GROUP}{$GENEID};
        $fasta_sequence_AAEGYPTI =~s/\|/\n/gis;

        if ( $MATCH_aegypti_albopictus_ortho_genes{$GENEID} ) {
            print TARGET_POPULATION "$fasta_sequence_AAEGYPTI\n"; 
        }


        ## >>> ADD ORTHOS <<<
        ## --------------------------------------------------------------------------------
        my $fasta_sequence_AALBOPICTUS = "";
        if ( $MATCH_aegypti_albopictus_ortho_genes{$GENEID} ) {
            my @fasta_line_AALBOPICTUS = split /\|/, $MATCH_aegypti_albopictus_ortho_genes{$GENEID};    ###/
               $fasta_sequence_AALBOPICTUS = $fasta_line_AALBOPICTUS[1];
               $fasta_sequence_AALBOPICTUS =~s/\t/\n/gis;

            print TARGET_POPULATION "$fasta_sequence_AALBOPICTUS\n";        
        }
        else {
            # seqs w/no orthologs
            #print OUTPUT_AFRICA_FULL "$aalbo_fasta_sequence\n";
        }

        close TARGET_POPULATION;
    }
}


print "All jobs done!\n\nCheck outfiles @ DIR:
/path/to/results/selection_dnds/aedes_aegypti_local_adaptation.yn_count_model/$POP_TARGET \n\n";





