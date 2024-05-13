#!/usr/bin/perl


use strict;
use warnings;






my %POPULATION_RAISD_GENES_counts_outlier;
my %POPULATION_RAISD_GENES_counts_genomic_feature;
my %POPULATION_RAISD_GENES_counts_strands;
my %POPULATION_GEO_INFORMATION;


my $population_information = "/scr/k80san3/ilozada/aedes_aegypti/NEE_paper/programs/scripts/running_on_slurm/aedes-aegypti.aaegL5.main.all_populations.general_information.byPopulation.updated_2023.tbl";

open POP_INFO, $population_information or die "CANNOT open the infile `POP_INFO` ~ $! \n";
my @population_table_info = ();
   @population_table_info = <POP_INFO>;
close POP_INFO;


foreach my $line (@population_table_info) {
    chomp $line;

    if ($line!~/^\#.+/ && $line!~/^Count.+/) {
        my @array     = split /\t/, $line;        
        my $continent = $array[2];
        my $popLine   = $array[$#array]; 
        
        my $group = "";
        if ($continent=~/Africa/) {  $group = "Africa";  }  else {  $group = "Out_of_Africa";  }

        my ($country,$population) = split /\./, $popLine;

        $POPULATION_GEO_INFORMATION{$country} = $group;
    }
}




# read directories
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my @main_directory = qw(/scr/k80san3/ilozada/aedes_aegypti/NEE_paper/programs/scripts/running_on_slurm/paper_NEE_rebuttal/raisd.changes_rebuttal.2024/raisd.redo_local.scripts_n_outputs/outputs/);

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
            
            if ($j=~/^selective_sweeps\.raisd\.(.+)\.(.+)\.(.+)\.fltrd_dataset\.outliers\.(.+)\.sweeps_to_refgenome\.vdate\.bed$/) {
		        my ($country,$population,$chromosome,$cutoff);
                $country    = $1;
                $population = $2;
                $chromosome = $3;
                $cutoff     = $4;


                open INFILE, $j or die "CANNOT open the file: `INFILE` $j - $!\n";

                # MAIN
                while (<INFILE>) {
                    chomp;
                    my $line = $_;

                    my ($chrm,$genomic_location,$genomic_location_cp,$startW,$endW,$VAR,$SFS,$LD,$mu_stat,$outlier_flag,$GFF_chrm,$GFF_start,$GFF_end,$GFF_geneID,$dot1,$GFF_strand,$GFF_db,$GFF_genomic_region,$dot2,$description) = split /\t/, $line;
                    my $outlier_info = "$chrm,$genomic_location,$startW,$endW,$VAR,$SFS,$LD,$mu_stat";

                    push( @{ $POPULATION_RAISD_GENES_counts_outlier{$cutoff}{$country}{$population}{$GFF_geneID} }, $outlier_info );
                    push( @{ $POPULATION_RAISD_GENES_counts_strands{$cutoff}{$country}{$population}{$GFF_geneID} }, $GFF_strand );
                    push( @{ $POPULATION_RAISD_GENES_counts_genomic_feature{$cutoff}{$country}{$population}{$GFF_geneID} }, $GFF_genomic_region );
                }

            }
        }
    }
}







my %hash_gene_counts_cutoff_05;
my %hash_gene_counts_cutoff_1;
my %hash_gene_counts_cutoff_05_intersections;
my %hash_gene_counts_cutoff_1_intersections;


my %genes_cutoff_05_OUT_of_AFRICA;
my %genes_cutoff_05_AFRICA;

my %genes_cutoff_1_OUT_of_AFRICA;
my %genes_cutoff_1_AFRICA;

my %genes_cutoff_05_OUT_of_AFRICA_outliers;
my %genes_cutoff_05_AFRICA_outliers;

my %genes_cutoff_1_OUT_of_AFRICA_outliers;
my %genes_cutoff_1_AFRICA_outliers;

my @full_geneset_list;
my @full_geneset_list_05;
my @full_geneset_list_1;


# OUTFILES SET 1: General info of each single gene outlier per population:
# ---------------------------------------------------------------------------------------------------------
my $result_sweeps_raisd_populations_thres_05_outliers = "/scr/k80san3/ilozada/aedes_aegypti/NEE_paper/programs/scripts/running_on_slurm/paper_NEE_rebuttal/raisd.changes_rebuttal.2024/raisd.redo_local.scripts_n_outputs/outputs/aedes_aegypti.populations.chrms123.selective_sweeps.raisd.parsed_results.thres_05.main.genes_n_outliers_info.v2.txt";
my $result_sweeps_raisd_populations_thres_1_outliers  = "/scr/k80san3/ilozada/aedes_aegypti/NEE_paper/programs/scripts/running_on_slurm/paper_NEE_rebuttal/raisd.changes_rebuttal.2024/raisd.redo_local.scripts_n_outputs/outputs/aedes_aegypti.populations.chrms123.selective_sweeps.raisd.parsed_results.thres_1.main.genes_n_outliers_info.v2.txt";

system "rm $result_sweeps_raisd_populations_thres_05_outliers  $result_sweeps_raisd_populations_thres_1_outliers";

open OUTFILE_GENES_N_OUTLIERS_ALL_05, ">>$result_sweeps_raisd_populations_thres_05_outliers" or die "CANNOT open the outfile `OUTFILE_GENES_N_OUTLIERS_ALL_05` ~ $!\n";
open OUTFILE_GENES_N_OUTLIERS_ALL_1, ">>$result_sweeps_raisd_populations_thres_1_outliers" or die "CANNOT open the outfile `OUTFILE_GENES_N_OUTLIERS_ALL_1` ~ $!\n";

print OUTFILE_GENES_N_OUTLIERS_ALL_05 "Group\tCutoff\tCountry\tPopulation\tGeneID\tTotal_outliers\tTotal_strands\tTotal_genome_features\tGenome_features\tStrands\tOutliers_info_merged\n";
print OUTFILE_GENES_N_OUTLIERS_ALL_1 "Group\tCutoff\tCountry\tPopulation\tGeneID\tTotal_outliers\tTotal_strands\tTotal_genome_features\tGenome_features\tStrands\tOutliers_info_merged\n";



# OUTFILES SET 2: Gene (outliers) distribution across populations as either shared or specific:
# ---------------------------------------------------------------------------------------------------------
my $result_sweeps_raisd_populations_thres_05_genes_PaP = "/scr/k80san3/ilozada/aedes_aegypti/NEE_paper/programs/scripts/running_on_slurm/paper_NEE_rebuttal/raisd.changes_rebuttal.2024/raisd.redo_local.scripts_n_outputs/outputs/aedes_aegypti.populations.chrms123.selective_sweeps.raisd.parsed_results.thres_05.main.genes_PaP_info.v2.txt";
my $result_sweeps_raisd_populations_thres_1_genes_PaP  = "/scr/k80san3/ilozada/aedes_aegypti/NEE_paper/programs/scripts/running_on_slurm/paper_NEE_rebuttal/raisd.changes_rebuttal.2024/raisd.redo_local.scripts_n_outputs/outputs/aedes_aegypti.populations.chrms123.selective_sweeps.raisd.parsed_results.thres_1.main.genes_PaP_info.v2.txt";

system "rm $result_sweeps_raisd_populations_thres_05_genes_PaP  $result_sweeps_raisd_populations_thres_1_genes_PaP";

open OUTFILE_GENES_PAP_05, ">>$result_sweeps_raisd_populations_thres_05_genes_PaP" or die "CANNOT open the outfile `OUTFILE_GENES_PAP_05` ~ $!\n";
open OUTFILE_GENES_PAP_1, ">>$result_sweeps_raisd_populations_thres_1_genes_PaP" or die "CANNOT open the outfile `OUTFILE_GENES_PAP_1` ~ $!\n";

print OUTFILE_GENES_PAP_05 "Cutoff\tGroup\tGeneID\tPopulation_counts_(Africa:Out_of_Africa)\tAfrica_populations\tOut_of_Africa_populations\n";
print OUTFILE_GENES_PAP_1 "Cutoff\tGroup\tGeneID\tPopulation_counts_(Africa:Out_of_Africa)\tAfrica_populations\tOut_of_Africa_populations\n";




for my $cutoff (keys %POPULATION_RAISD_GENES_counts_outlier) {
    for my $country (keys %{ $POPULATION_RAISD_GENES_counts_outlier{$cutoff} } ) {

        # population group:
        my $pop_group = "";
        if ( $POPULATION_GEO_INFORMATION{$country} ) {  $pop_group = $POPULATION_GEO_INFORMATION{$country};  } 

        for my $population (keys %{ $POPULATION_RAISD_GENES_counts_outlier{$cutoff}{$country} } ) {

            my @genes_per_population_thres_05;
            my @genes_per_population_thres_1;
            
            for my $geneID (keys %{ $POPULATION_RAISD_GENES_counts_outlier{$cutoff}{$country}{$population} } ) {

                # Parse data to make FULL REPORT:
                # ---------------------------------------------------------------------------------------------------------------------------------------
                my @outliers_unique = uniq( @{ $POPULATION_RAISD_GENES_counts_outlier{$cutoff}{$country}{$population}{$geneID} } );
                my $outliers_merged = join "|", @outliers_unique;
                my $total_outliers = scalar (@outliers_unique);

                my @strand_unique = "";  my $strand_merged = "";  my $total_strands = "";

                # section: strands
                if ( $POPULATION_RAISD_GENES_counts_strands{$cutoff}{$country}{$population}{$geneID} ) {
                    @strand_unique = uniq( @{ $POPULATION_RAISD_GENES_counts_strands{$cutoff}{$country}{$population}{$geneID} } );
                    $strand_merged = join "|", @strand_unique;
                    $total_strands = scalar (@strand_unique);
                }


                my @genomfeats_unique = "";  my $genomfeats_merged = "";  my $total_genomfeats  = "";

                # section: genomic features
                if ( $POPULATION_RAISD_GENES_counts_genomic_feature{$cutoff}{$country}{$population}{$geneID} ) {
                    @genomfeats_unique = uniq( @{ $POPULATION_RAISD_GENES_counts_genomic_feature{$cutoff}{$country}{$population}{$geneID} } );
                    $genomfeats_merged = join "|", @genomfeats_unique;
                    $total_genomfeats = scalar (@genomfeats_unique);
                }


                # MAIN: CUTOFF 0.05
                if ($cutoff=~/thres_05/) {
                    print OUTFILE_GENES_N_OUTLIERS_ALL_05 "$pop_group\t0.05\t$country\t$population\t$geneID\t$total_outliers\t$total_strands\t$total_genomfeats\t$genomfeats_merged\t$strand_merged\t$outliers_merged\n";
                    push(@genes_per_population_thres_05,$geneID);
                    
                }
                elsif ($cutoff=~/thres_1/) {
                    print OUTFILE_GENES_N_OUTLIERS_ALL_1 "$pop_group\t1.0\t$country\t$population\t$geneID\t$total_outliers\t$total_strands\t$total_genomfeats\t$genomfeats_merged\t$strand_merged\t$outliers_merged\n";
                    push(@hash_gene_counts_cutoff_1,$geneID);
                    
                }
                else {
                    #print "$cutoff\t$country\t$population\t$geneID\n";
                }
                # ---------------------------------------------------------------------------------------------------------------------------------------




                # ---------------------------------------------------------------------------------------------------------------------------------------
                # Get genes for Africa and Out of Africa separately:
                # ---------------------------------------------------------------------------------------------------------------------------------------
                if ($pop_group=~/^Africa$/) {
                    
                    if ($cutoff=~/thres_05/)   {  push( @{ $genes_cutoff_05_AFRICA{$geneID} }, $population );  push(@full_geneset_list_05,$geneID);  push( @{ $genes_cutoff_05_AFRICA_outliers{$geneID} }, $outliers_merged );  }
                    elsif ($cutoff=~/thres_1/) {  push( @{ $genes_cutoff_1_AFRICA{$geneID} }, $population );   push(@full_geneset_list_1,$geneID);   push( @{ $genes_cutoff_1_AFRICA_outliers{$geneID} }, $outliers_merged );  }
                    else { print "WARNING_AFRICA!!!$geneID\t$population\n";}
                }
                elsif ($pop_group=~/^Out_of_Africa$/) {
                    
                    if ($cutoff=~/thres_05/)   { push( @{ $genes_cutoff_05_OUT_of_AFRICA{$geneID} }, $population ); push(@full_geneset_list_05,$geneID);  push( @{ $genes_cutoff_05_OUT_of_AFRICA_outliers{$geneID} }, $outliers_merged );  }
                    elsif ($cutoff=~/thres_1/) { push( @{ $genes_cutoff_1_OUT_of_AFRICA{$geneID} }, $population );  push(@full_geneset_list_1,$geneID);   push( @{ $genes_cutoff_1_OUT_of_AFRICA_outliers{$geneID} }, $outliers_merged );  }
                    else { print "WARNING_OUT_OF_AFRICA!!!$geneID\t$population\n";}
                }
                else { print "$pop_group\t$country\t$population\t$total_genes_05\n"; }
                # ---------------------------------------------------------------------------------------------------------------------------------------
                
                push(@full_geneset_list,$geneID);
            }

            my $total_genes_05 = scalar( uniq(@genes_per_population_thres_05) );
            my $total_genes_1  = scalar( uniq(@genes_per_population_thres_1) );

            if ($cutoff=~/thres_05/) {
                $hash_gene_counts_cutoff_05{$population} = $total_genes_05;
            }
            if ($cutoff=~/thres_1/) {
                $hash_gene_counts_cutoff_1{$population} = $total_genes_1;
            }

        }
    }
}
close OUTFILE_GENES_N_OUTLIERS_ALL_05;
close OUTFILE_GENES_N_OUTLIERS_ALL_1;





# --------------------------------------------------------------------------------------------------------------------------------------------
# shared genes at the 0.05% cutoff
# --------------------------------------------------------------------------------------------------------------------------------------------

my @full_geneset_list_05_unique = uniq(@full_geneset_list_05);

foreach my $geneID (@full_geneset_list_05_unique) {
    chomp $geneID;
    #print "$geneID\n"; # total genes 6927

    if ( $genes_cutoff_05_AFRICA{$geneID} && $genes_cutoff_05_OUT_of_AFRICA{$geneID} && $genes_cutoff_05_AFRICA_outliers{$geneID} && $genes_cutoff_05_OUT_of_AFRICA_outliers{$geneID} ) {
        # SHARED: 2658
        my $PopList_OoA   = join "|", uniq( @{ $genes_cutoff_05_OUT_of_AFRICA{$geneID} } );
        my $PopList_AFR   = join "|", uniq( @{ $genes_cutoff_05_AFRICA{$geneID} } );
        my $PopCounts_OoA = scalar (uniq( @{ $genes_cutoff_05_OUT_of_AFRICA{$geneID} } ));
        my $PopCounts_AFR = scalar (uniq( @{ $genes_cutoff_05_AFRICA{$geneID} } ));

        my $OutliersList_OoA   = join "|", uniq( @{ $genes_cutoff_05_OUT_of_AFRICA_outliers{$geneID} } );
        my $OutliersList_AFR   = join "|", uniq( @{ $genes_cutoff_05_AFRICA_outliers{$geneID} } );
        my $OutliersList_BOTH_MERGED = "[OoA\:$OutliersList_OoA]\[AFR\:$OutliersList_AFR]";

        print OUTFILE_GENES_PAP_05 "0.05\%\tShared\t$geneID\t$PopCounts_AFR\:$PopCounts_OoA\t$PopList_AFR\t$PopList_OoA\t$OutliersList_BOTH_MERGED\n";
    }
    elsif ( $genes_cutoff_05_AFRICA{$geneID} && !$genes_cutoff_05_OUT_of_AFRICA{$geneID} && $genes_cutoff_05_AFRICA_outliers{$geneID} && !$genes_cutoff_05_OUT_of_AFRICA_outliers{$geneID} ) {
        # AFR: 3375
        my $PopList_OoA   = "na";
        my $PopList_AFR   = join "|", uniq( @{ $genes_cutoff_05_AFRICA{$geneID} } );
        my $PopCounts_OoA = "na";
        my $PopCounts_AFR = scalar (uniq( @{ $genes_cutoff_05_AFRICA{$geneID} } ));

        my $OutliersList_AFR   = join "|", uniq( @{ $genes_cutoff_05_AFRICA_outliers{$geneID} } );

        print OUTFILE_GENES_PAP_05 "0.05\%\tAfrica_specific\t$geneID\t$PopCounts_AFR\:$PopCounts_OoA\t$PopList_AFR\t$PopList_OoA\t$OutliersList_AFR\n";
    }
    elsif ( !$genes_cutoff_05_AFRICA{$geneID} && $genes_cutoff_05_OUT_of_AFRICA{$geneID} && !$genes_cutoff_05_AFRICA_outliers{$geneID} && $genes_cutoff_05_OUT_of_AFRICA_outliers{$geneID} ) {
        # OoA: 894
        my $PopList_OoA   = join "|", uniq( @{ $genes_cutoff_05_OUT_of_AFRICA{$geneID} } );
        my $PopList_AFR   = "na";
        my $PopCounts_OoA = scalar (uniq( @{ $genes_cutoff_05_OUT_of_AFRICA{$geneID} } ));
        my $PopCounts_AFR = "na";

        my $OutliersList_OoA   = join "|", uniq( @{ $genes_cutoff_05_OUT_of_AFRICA_outliers{$geneID} } );

        print OUTFILE_GENES_PAP_05 "0.05\%\tOut_of_Africa_specific\t$geneID\t$PopCounts_AFR\:$PopCounts_OoA\t$PopList_AFR\t$PopList_OoA\t$OutliersList_OoA\n";
    }
    else {
        print "WARNINIG TOP 0.05%!!! $geneID\t|NA|\t|NA|\n";
    }
}
close OUTFILE_GENES_PAP_05;





# --------------------------------------------------------------------------------------------------------------------------------------------
# shared genes at the 1.0% cutoff
# --------------------------------------------------------------------------------------------------------------------------------------------

my @full_geneset_list_1_unique = uniq(@full_geneset_list_1);

foreach my $geneID (@full_geneset_list_1_unique) {
    chomp $geneID;

    if ( $genes_cutoff_1_AFRICA{$geneID} && $genes_cutoff_1_OUT_of_AFRICA{$geneID} && $genes_cutoff_1_AFRICA_outliers{$geneID} && $genes_cutoff_1_OUT_of_AFRICA_outliers{$geneID} ) {
        my $PopList_OoA   = join "|", uniq( @{ $genes_cutoff_1_OUT_of_AFRICA{$geneID} } );
        my $PopList_AFR   = join "|", uniq( @{ $genes_cutoff_1_AFRICA{$geneID} } );
        my $PopCounts_OoA = scalar (uniq( @{ $genes_cutoff_1_OUT_of_AFRICA{$geneID} } ));
        my $PopCounts_AFR = scalar (uniq( @{ $genes_cutoff_1_AFRICA{$geneID} } ));

        my $OutliersList_OoA   = join "|", uniq( @{ $genes_cutoff_1_OUT_of_AFRICA_outliers{$geneID} } );
        my $OutliersList_AFR   = join "|", uniq( @{ $genes_cutoff_1_AFRICA_outliers{$geneID} } );
        my $OutliersList_BOTH_MERGED = "[OoA\:$OutliersList_OoA]\[AFR\:$OutliersList_AFR]";

        print OUTFILE_GENES_PAP_1 "1.0\%\tShared\t$geneID\t$PopCounts_AFR\:$PopCounts_OoA\t$PopList_AFR\t$PopList_OoA\t$OutliersList_BOTH_MERGED\n";
    }
    elsif ( $genes_cutoff_1_AFRICA{$geneID} && !$genes_cutoff_1_OUT_of_AFRICA{$geneID} && $genes_cutoff_1_AFRICA_outliers{$geneID} && !$genes_cutoff_1_OUT_of_AFRICA_outliers{$geneID} ) {
        my $PopList_OoA   = "na";
        my $PopList_AFR   = join "|", uniq( @{ $genes_cutoff_1_AFRICA{$geneID} } );
        my $PopCounts_OoA = "na";
        my $PopCounts_AFR = scalar (uniq( @{ $genes_cutoff_1_AFRICA{$geneID} } ));

        my $OutliersList_AFR   = join "|", uniq( @{ $genes_cutoff_1_AFRICA_outliers{$geneID} } );

        print OUTFILE_GENES_PAP_1 "1.0\%\tAfrica_specific\t$geneID\t$PopCounts_AFR\:$PopCounts_OoA\t$PopList_AFR\t$PopList_OoA\t$OutliersList_AFR\n";
    }
    elsif ( !$genes_cutoff_1_AFRICA{$geneID} && $genes_cutoff_1_OUT_of_AFRICA{$geneID} && !$genes_cutoff_1_AFRICA_outliers{$geneID} && $genes_cutoff_1_OUT_of_AFRICA_outliers{$geneID} ) {
        my $PopList_OoA   = join "|", uniq( @{ $genes_cutoff_1_OUT_of_AFRICA{$geneID} } );
        my $PopList_AFR   = "na";
        my $PopCounts_OoA = scalar (uniq( @{ $genes_cutoff_1_OUT_of_AFRICA{$geneID} } ));
        my $PopCounts_AFR = "na";

        my $OutliersList_OoA   = join "|", uniq( @{ $genes_cutoff_1_OUT_of_AFRICA_outliers{$geneID} } );

        print OUTFILE_GENES_PAP_1 "1.0\%\tOut_of_Africa_specific\t$geneID\t$PopCounts_AFR\:$PopCounts_OoA\t$PopList_AFR\t$PopList_OoA\t$OutliersList_OoA\n";
    }
    else {
        print "WARNINIG TOP 1%!!! $geneID\t|NA|\t|NA|\n";
    }
}
close OUTFILE_GENES_PAP_1;











sub uniq {
  my %seen;
  return grep { !$seen{$_}++ } @_;
}

