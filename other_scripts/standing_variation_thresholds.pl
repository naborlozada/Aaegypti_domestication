#!/usr/bin/perl

# ------------------------------------------------ #
#      *** Alejandro Nabor Lozada Chavez ***       #
# Questions: alejandro.chavez@unipv.it             #
#            nabor.lozada@gmail.com                #
# ------------------------------------------------ #

use strict;
use warnings;
use IO::Zlib;
use Compress::Zlib;
use List::MoreUtils qw(uniq);





# Variables
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my $localdate = localtime();

my %pnps_byGene_byPop;
my %mine_Pops_list;
my @gene_ids_allPushed;

my %Standing_Variation_genes_extended; 

my $pnps_whole_genome="aaeg.pops.selection.pnps_info.whole_genome.txt";

open PNPS_WHOLE_GENOME, "$pnps_whole_genome" or die "CANNOT open the INFILE of pNpS PNPS_WHOLE_GENOME $pnps_whole_genome ~ $!\n";

my @pnps_whole_genome_array = ();
   @pnps_whole_genome_array = <PNPS_WHOLE_GENOME>;
close PNPS_WHOLE_GENOME;

foreach my $line (@pnps_whole_genome_array) {
    chomp $line;
    if ($line!~/^#.+/ && $line!~/^Country/) {
        my @split_info    = split /\t/, $line; 
        my $population    = $split_info[1];
        my $geneID        = $split_info[3];
        my $pnps_score    = $split_info[16];
        my $selectionType = $split_info[20];

        push(@gene_ids_allPushed,$geneID);

        # define selection in numeric code: from R script
        my $neutrality1_rlx = 0.80;   my $neutrality2_rlx = 1.20;   

        # standing variation range:  0.80 -- 1.20
        if ($pnps_score!~/^NA$/ && $pnps_score!~/^Inf$/ && $pnps_score >= $neutrality1_rlx && $pnps_score <= $neutrality2_rlx) {
            #print "$geneID\t$pnps_score\t$selectionMode\t$population\n";
            push( @{ $Standing_Variation_genes_extended{$geneID} }, $population );
        }
    }
}





my @unique_geneID = uniq(@gene_ids_allPushed);
my $total_unique_geneID = scalar @unique_geneID;
#print "[$total_unique_geneID]\n";


# Populate an array with our populations
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my @defined_population_list = ("thies","ngoye","bantata","kedougou","mindin","pk10","dori","ouagadougou","ouahigouya","larabanga","boabengFiema","kintampo","kumasi","awka","benoue","franceville","libreville","lopeVillage","bundibugyo","entebbe","karenga","kichwamba","virhembe","kakamega","Mbarakani_village","ganda","arabuko","kayaBomu","kwale","shimbaHills","rabai_selv","rabai_dom","tapachula","bebeouro","santarem","jeddah_al-rawabi","samut_sakhon","bangkok","tafuna_village","zac_panda");
my %population_information;

my $infile_info = "/home/tigerpv_das/users_data/Nabor/Aedes-aegypti/MAIN/analysis/selection/kaks/matrices_k2p_paml/infile_popnames_to_match.txt";
open INFILE_INFO, "$infile_info" or die "CANNOT open the infile $infile_info ~ $!\n";
my @popNames = ();
   @popNames = <INFILE_INFO>;
close INFILE_INFO;

foreach my $line (@popNames){
    chomp $line;
    #print "$line\n";
    my ($nr,$region,$population_lower,$population,$code) = split /\t/, $line;
    $population_information{$population_lower} = $region;
}



# Make outfiles (Supplementary Table 18)
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my $standing_variation_geneset_25prcnt = "aedes_aegypti.stading_variation.whole_genome.25prcnt_threshold.geneset_list.txt";
my $standing_variation_geneset_50prcnt = "aedes_aegypti.stading_variation.whole_genome.50prcnt_threshold.geneset_list.txt";
my $standing_variation_geneset_75prcnt = "aedes_aegypti.stading_variation.whole_genome.75prcnt_threshold.geneset_list.txt";

open STADVAR_25PERCENT, ">>$standing_variation_geneset_25prcnt" or die "CANNOT open the OUTFILE $standing_variation_geneset_25prcnt ~ $!\n";
open STADVAR_50PERCENT, ">>$standing_variation_geneset_50prcnt" or die "CANNOT open the OUTFILE $standing_variation_geneset_50prcnt ~ $!\n";
open STADVAR_75PERCENT, ">>$standing_variation_geneset_75prcnt" or die "CANNOT open the OUTFILE $standing_variation_geneset_75prcnt ~ $!\n";

print STADVAR_25PERCENT "# version: $localdate\n";
print STADVAR_50PERCENT "# version: $localdate\n";
print STADVAR_75PERCENT "# version: $localdate\n";

print STADVAR_25PERCENT "# This file contains a list of genes under STADING VARIATION based on a relaxed threshold of pN/pS ratio estimated using K2P approach\n";
print STADVAR_50PERCENT "# This file contains a list of genes under STADING VARIATION based on a relaxed threshold of pN/pS ratio estimated using K2P approach\n";
print STADVAR_75PERCENT "# This file contains a list of genes under STADING VARIATION based on a relaxed threshold of pN/pS ratio estimated using K2P approach\n";

print STADVAR_25PERCENT "Percentage\tGene_ID\tTotal_Nr_Populations\tNr_African_populations\tNr_out-of-Africa_populations\tAll_populations\n";
print STADVAR_50PERCENT "Percentage\tGene_ID\tTotal_Nr_Populations\tNr_African_populations\tNr_out-of-Africa_populations\tAll_populations\n";
print STADVAR_75PERCENT "Percentage\tGene_ID\tTotal_Nr_Populations\tNr_African_populations\tNr_out-of-Africa_populations\tAll_populations\n";


# push( @{ $Standing_Variation_genes_extended{$geneID} }, $population );
my $total_count_StandingVariation_genes_25PopPrcnt = 0;
my $total_count_StandingVariation_genes_50PopPrcnt = 0;  
my $total_count_StandingVariation_genes_75PopPrcnt = 0;
my $total_genes_falling_in_Stading_variation = 0;

for my $genes_StandVar (keys %Standing_Variation_genes_extended) {
    
    my $total_Nr_populations = scalar @{ $Standing_Variation_genes_extended{$genes_StandVar} };
    my $all_populations      = join " ", @{ $Standing_Variation_genes_extended{$genes_StandVar} };

    my ($region,$popCount_Africa,$popCount_OutOfAfrica) = ("",0,0);
    
    foreach my $population ( @{ $Standing_Variation_genes_extended{$genes_StandVar} } ) {
        chomp $population;
        if ( $population_information{$population} ) {   $region = $population_information{$population};   }
        
        if ($region=~/^africa$/)           {  $popCount_Africa++;  }
        elsif ($region=~/^out_of_Africa$/) {  $popCount_OutOfAfrica++;  }
        else {  print "WARNING: $population \n";  }
    }
    
    #print "$genes_StandVar\t$total_Nr_populations\t[$popCount_Africa *** $popCount_OutOfAfrica]\t$all_populations\n";
    $total_genes_falling_in_Stading_variation++;
    
    # Set of genes under standing variation in at least 25% of all populations
    if ($total_Nr_populations>=10) {
        $total_count_StandingVariation_genes_25PopPrcnt++;
        print STADVAR_25PERCENT "25\t$genes_StandVar\t$total_Nr_populations\t$popCount_Africa\t$popCount_OutOfAfrica\t$all_populations\n";
    }
    # Set of genes under standing variation in at least 50% of all populations
    if ($total_Nr_populations>=20) {
        $total_count_StandingVariation_genes_50PopPrcnt++;
        print STADVAR_50PERCENT "50\t$genes_StandVar\t$total_Nr_populations\t$popCount_Africa\t$popCount_OutOfAfrica\t$all_populations\n";
    }
    # Set of genes under standing variation in at least 75% of all populations
    if ($total_Nr_populations>=30) {
        $total_count_StandingVariation_genes_75PopPrcnt++;
        print STADVAR_75PERCENT "75\t$genes_StandVar\t$total_Nr_populations\t$popCount_Africa\t$popCount_OutOfAfrica\t$all_populations\n";
    }

}


# additional annotation for 25%
# ---------------------------------------------------------------------------------------------------
# percentage of standing variation in the genome
my $percentage_whole_genome_25prcnt = ($total_count_StandingVariation_genes_25PopPrcnt * 100) / 14677; 

# percentage of genes that only fall in standing variation threshols 
my $percentage_standVar_genes_25prcnt = ($total_count_StandingVariation_genes_25PopPrcnt * 100) / $total_genes_falling_in_Stading_variation; 

print STADVAR_25PERCENT "Population minimum at least: 10 populations (25%)\n";
print STADVAR_25PERCENT "Normalized from total genome [$total_count_StandingVariation_genes_25PopPrcnt | 14677]:";
print STADVAR_25PERCENT "$total_count_StandingVariation_genes_25PopPrcnt ---> $percentage_whole_genome_25prcnt %\n";

print STADVAR_25PERCENT "Normalized based on the total set of genes in Standing Variation [$total_count_StandingVariation_genes_25PopPrcnt | $total_genes_falling_in_Stading_variation]:";
print STADVAR_25PERCENT "$total_count_StandingVariation_genes_25PopPrcnt ---> $percentage_standVar_genes_25prcnt % \n\n";




# additional annotation for 50%
# ---------------------------------------------------------------------------------------------------
# percentage of standing variation in the genome
my $percentage_whole_genome_50prcnt = ($total_count_StandingVariation_genes_50PopPrcnt * 100) / 14677; 

# percentage of genes that only fall in standing variation threshols 
my $percentage_standVar_genes_50prcnt = ($total_count_StandingVariation_genes_50PopPrcnt * 100) / $total_genes_falling_in_Stading_variation; 

print STADVAR_50PERCENT "Population minimum at least: 20 populations (50%)\n";
print STADVAR_50PERCENT "Normalized from total genome [$total_count_StandingVariation_genes_50PopPrcnt | 14677]:";
print STADVAR_50PERCENT "$total_count_StandingVariation_genes_50PopPrcnt ---> $percentage_whole_genome_50prcnt %\n";

print STADVAR_50PERCENT "Normalized based on the total set of genes in Standing Variation [$total_count_StandingVariation_genes_50PopPrcnt | $total_genes_falling_in_Stading_variation]:";
print STADVAR_50PERCENT "$total_count_StandingVariation_genes_50PopPrcnt ---> $percentage_standVar_genes_50prcnt % \n\n";



# additional annotation for 75%
# ---------------------------------------------------------------------------------------------------
# percentage of standing variation in the genome
my $percentage_whole_genome_75prcnt = ($total_count_StandingVariation_genes_75PopPrcnt * 100) / 14677; 

# percentage of genes that only fall in standing variation threshols 
my $percentage_standVar_genes_75prcnt = ($total_count_StandingVariation_genes_75PopPrcnt * 100) / $total_genes_falling_in_Stading_variation; 

print STADVAR_75PERCENT "Population minimum at least: 30 populations (75%)\n";
print STADVAR_75PERCENT "Normalized from total genome [$total_count_StandingVariation_genes_75PopPrcnt | 14677]:";
print STADVAR_75PERCENT "$total_count_StandingVariation_genes_75PopPrcnt ---> $percentage_whole_genome_75prcnt %\n";

print STADVAR_75PERCENT "Normalized based on the total set of genes in Standing Variation [$total_count_StandingVariation_genes_75PopPrcnt | $total_genes_falling_in_Stading_variation]:";
print STADVAR_75PERCENT "$total_count_StandingVariation_genes_75PopPrcnt ---> $percentage_standVar_genes_75prcnt % \n\n";


close STADVAR_25PERCENT;
close STADVAR_50PERCENT;
close STADVAR_75PERCENT;


# END
