#!/usr/bin/perl

# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada Chavez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #

## Description:
## ===========================================================================================================================================================================================
## This script gets all the statistics of D,S,Dn/Ds,Ds,Dn calculated with CODEML from the PAML package and performs three main tasks:
##      1. Creates an output file with filtered calculations (D,S,Dn/Ds,Ds,Dn) for each gene and rename it with the source population: name of country, population, chromosome, and gene ID.
##      2. For each gene, basic descriptive statistics are performed (mean, stdev, max, min, and mode) 
##      3. Save all information in an output file.
## ===========================================================================================================================================================================================


use strict;
use warnings;
use IO::Zlib;
use Compress::Zlib;
use List::MoreUtils qw(uniq);
use Statistics::Basic qw(:all);
use Statistics::Descriptive;
#use diagnostics;




# Variables
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my $localdate      = localtime();
my $population_DIR = "";    
my $main_DIR       = "";
my $main_FINAL_DIR = "";

my %HASH_CODEML_DNDS_SCORES;
my %HASH_CODEML_PAIR_NAMES;




# Make outfile
# ---------------------------------------------------------------------------------------------------------------------------------------------------
# aaeg.pops.selection.kaks.main_table.kaks_info.reference_gene_families.txt
my $outfile_codeml_dnds_scores="aedes_aegypti.selection.dnds.codeml_stats.whole_genome.all_populations_genes.output.txt";
system "rm $outfile_codeml_dnds_scores";
open OUTFILE_CODEML_DNDS_SCORES, ">>$outfile_codeml_dnds_scores" or die "CANNOT open the OUTFILE $outfile_codeml_dnds_scores ~ $! \n";
print OUTFILE_CODEML_DNDS_SCORES "# version: $localdate\n";
print OUTFILE_CODEML_DNDS_SCORES "# Info: Basic descriptive statistics for pairwise comparisons between sequences of a single population. Each single gene alignment was obtained only for those protein coding genes from a VCF file that present at least 1 SNP, and each population might have a different set of genes.\n";  
print OUTFILE_CODEML_DNDS_SCORES "# Protocol: [a] SNP sequences were transformed into fasta format using the 'vcf2fasta' script (REF). [b] Sequences with ambiguous nucleotides and/or with strange characters (ie. /?) were removed (BASH scripting). [c] Codon alignments were created (pal2nal.pl) after translate the original fasta sequences (transeq, EMBOSS). [d] Ratio of Dn/Ds was estimated (codeml, PAML) defining individial control files (BASH scripting), which all of them were submmited via 'parallel' GNU function (EOS + NIRV2 cluster nodes).\n";

print OUTFILE_CODEML_DNDS_SCORES "gene_file\tFILE_STATUS\tdnds_mean\tdnds_stdev\tdnds_min\tdnds_max\tdnds_mode\tCOUNTRY\tPOPULATION\tCHROMOSOME\toutfile\n";
# ---------------------------------------------------------------------------------------------------------------------------------------------------





# read directories
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my @main_directory = qw(/home/nlozada/aedes_aegpypti/analyses/dnds/results_from_eos/);

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

            if ($j!~/^\.$/ && $j!~ /^\.\.$/ && $j!~/^\.listing$/ && $j!~/readme/ && $j!~/.+\.txt$/ && $j=~/^(.+)\.(.+)\.(.+)\_CDS/) {
		        my ($country,$population,$chromosome);
                $country    = $1;
                $population = $2;
                $chromosome = $3;
                #print "$country\t$population\t$chromosome\t$j\n";
                my $subpwd = "/home/nlozada/aedes_aegpypti/analyses/dnds/results_from_eos/$j/";
		            $main_DIR=$j;
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

            if ($j !~ /^\.$/ && $j!~ /^\.\.$/ && $j!~ /^\.listing$/ && $j!~/.+\.fas$/ && $j!~/.+\.txt$/ && $j!~/.+\.tab$/ && $j!~/snap\.background/ && $j!~/snap\.summary/ && $j=~/^(.+)\.nts\.fixedN\.codons\.codeml\.ctldir$/) {
		        my $subpwd2 = "/home/nlozada/aedes_aegpypti/analyses/dnds/results_from_eos/$main_DIR/$j/";        ### aedes-aegypti.admixture_groups.geno02
                chdir $subpwd2;
                opendir SUBSUBDIR, $subpwd2 or die "CANNOT open the SUB-DIR $subpwd2 - $! \n";
		            $main_FINAL_DIR=$subpwd2;

                my @files2=();
		            @files2=readdir SUBSUBDIR;
                &GoGoDown3(@files2);
                closedir SUBSUBDIR;
            }
        }
    }
}


# >>> GO-GO-GO to NEXT SUBDIR >>>
sub GoGoDown3{

    my @files2= @_;
    if ($#files2 > 1){
        @files2 = @files2[0..$#files2];
        my $files = @files2;
        my @files_sorted = sort {lc($a) cmp lc($b)} @files2;

        foreach my $local_file (@files_sorted){
            chomp $local_file;

            my $population_chrom = "";
            if ($local_file=~/^(.+)\.nts\.fixedN\.codons\.out$/) {
                my $gene_file=$1;
                my @array2split   = split /\//, $main_FINAL_DIR;
                my $pop_full_info = $array2split[7];
                my $outfile       = "$pop_full_info\.$gene_file\.paml\.codeml\.output\.txt";
                #print "$main_FINAL_DIR\t$outfile\t$local_file\t[[$pop_full_info]]\n"; 
                
                my ($COUNTRY,$POPULATION,$CHROMOSOME) = split /\./, $pop_full_info;
                $CHROMOSOME=~ s/\_CDS$//gis;
                open CODEML_OUTPUT, $local_file or die "CANNOT open the file $local_file - $!\n";
                my $index_scores = 0;
                my $index_names  = 0;
                my $index_lines  = 0;
                my @dnds_gene_scores=();

                # MAIN
                # ------------------------------------------------------------------------------------------
                while (<CODEML_OUTPUT>) {
                    chomp;
                    my $line = $_;
                    $index_lines++;
                    
                    my ($S,$N,$DnDs,$dN,$dS) = (0,0,0,0,0);
                    my ($query,$subject) = ("","");
                    # searching structure
                    #   t= 0.0000  S=    41.1  N=   126.9  dN/dS=  1.3546  dN = 0.0000  dS = 0.0000
                    #   t= 0.0000  S=   254.0  N=  1051.0  dN/dS=  0.0010  dN = 0.0000  dS = 0.0000
                    #   t= 0.0000  S=   254.0  N=  1051.0  dN/dS=  0.0010  dN = 0.0000  dS = 0.0000

                    if ($line=~/^t\=/ && $line=~/dN\/dS\=/) {
                        $index_scores++;
                        my @splits = split /\s+/, $line;
                        $S    = $splits[3];
                        $N    = $splits[5];
                        $DnDs = $splits[7];
                        $dN   = $splits[10];
                        $dS   = $splits[13];
                        
                        #print "$index_scores\t[$S  $N  $DnDs  $dN  $dS]\t$line\n";
                        push(@dnds_gene_scores,$DnDs);                     
                    }
                }

                # Defina variables
                my ($dnds_mean,$dnds_stdev,$dnds_min,$dnds_max,$dnds_mode,$FILE_STATUS);

                # Do basic stastistics
                if ( $index_lines != 0 ) {
                    $FILE_STATUS = "COMPLETE";
                    #my @data = (1..10);
                    my $stat = Statistics::Descriptive::Full->new();
                       $stat->add_data(@dnds_gene_scores);
                       $dnds_mean  = $stat->mean();
                       $dnds_stdev = $stat->standard_deviation();
                       $dnds_min   = $stat->min();
                       $dnds_max   = $stat->max();
                       $dnds_mode  = $stat->mode();
                }
                else {
                    $FILE_STATUS = "EMPTY";
                    $dnds_mean="NA";  $dnds_stdev="NA";  $dnds_min="NA";  $dnds_max="NA";  $dnds_mode="NA";
                }

                print OUTFILE_CODEML_DNDS_SCORES "$gene_file\t$FILE_STATUS\t$dnds_mean\t$dnds_stdev\t$dnds_min\t$dnds_max\t$dnds_mode\t$COUNTRY\t$POPULATION\t$CHROMOSOME\t$outfile\n";
                # ------------------------------------------------------------------------------------------
            }
        }
    }
}
close OUTFILE_CODEML_DNDS_SCORES;



# END




