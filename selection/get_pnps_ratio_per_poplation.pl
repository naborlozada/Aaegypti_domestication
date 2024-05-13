#!/usr/bin/perl

use strict;
use warnings;



# Variables
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my $localdate      = localtime();
my $population_DIR = "";    
my $main_DIR       = "";
my $main_FULL_DIR  = "";

my %MAJOR_POPULATION_GROUPS;
my %aegypti_consensus_sequences;

my %HASH_ORTHOLOGS_GENES;




# prev: aedes_aegypti_local_adaptation.yn_count_model
my $main_OUTPUT_DIR = "/path/to/results/selection_dnds/parsed_results";
my $pal2nal = "/path/to/software_scripts/pal2nal.v14/pal2nal.pl";
my $fas2axt = "/path/to/results/selection_dnds/scripts/selection.GLOBAL.pairwise/parseFastaIntoAXT.pl";


my $outfile_dnds_selection_YN00_SUPPL_DATA_7 = "aedes_aegypti.populations.local_selection.dNdS.main_table.whole_genome.gene_scores.PAML_basic_info.SUPPL_DATA_7.txt";

system "rm $main_OUTPUT_DIR/$outfile_dnds_selection_YN00_SUPPL_DATA_7";

open OUTFILE_SUPPLEMENTARY_DATA_SEVEN, ">>$main_OUTPUT_DIR/$outfile_dnds_selection_YN00_SUPPL_DATA_7" or die "CANNOT open the OUTFILE_SUPPLEMENTARY_DATA_SEVEN $outfile_dnds_selection_YN00_SUPPL_DATA_7 ~ $! \n";
print OUTFILE_SUPPLEMENTARY_DATA_SEVEN "LOCATION\tCOUNTRY\tPOPULATION\tAeAEGYPTI_GENE_ID\tAeALBOPICTUS_GENE_ID\tSELECTION_TYPE\tNON_SYNONYMOUS_N\tSYNONYMOUS_D\tOMEGA\tNON_SYNONYMOUS_N_RATIO\tSYNONYMOUS_D_RATIO\tKAPPA\tt_VALUE\n";


  


my $ortholgues_genes_list = "/path/to/results/selection_dnds/parsed_results/aedes_aegypti.orthologues_genes.ae_aegypti__vs__ae_albopictus.whole_genome_comparison.two_assemblies.table.txt";

open ORTHO_GENES_LIST, "$ortholgues_genes_list" or die "CANNOT open the INFILE `ORTHO_GENES_LIST` ~ $!\n";
my @orthos_genes = ();
   @orthos_genes = <ORTHO_GENES_LIST>;
close ORTHO_GENES_LIST;


foreach my $gene_line (@orthos_genes) {
    chomp $gene_line;
    if ($gene_line!~/^#/ && $gene_line!~/^Aedes_aegypti_ID/) {
        my ($aegypti,$albopictus,$assembly) = split /\t/, $gene_line;
        $HASH_ORTHOLOGS_GENES{$aegypti} = $albopictus;
    }
}







# groups
# ---------------------------------------------------------------------------------------------------------------------------------------------------
my $groups_file = "path/to/results/selection_dnds/scripts/selection.GLOBAL.pairwise/aedes_aegypti.populations_information.groups_n_geography.tbl.txt";
open AAEL_GROUPS, "$groups_file" or die "CANNOT open the infile AAEL_GROUPS ~ $!\n";
my @groups = ();
   @groups = <AAEL_GROUPS>;
close AAEL_GROUPS;


foreach my $line (@groups) {
    chomp $line;

    if ( $line!~/^\#.+/ && $line!~/^Nr.+/ ) {
        my @split_info = split /\t/, $line;
        my $POPdir = $split_info[1];
        my $GROUP1 = $split_info[2];    ## FULL: AFRICA / OoA
        my $GROUP2 = $split_info[3];    ## FULL REGIONS: WEST / EAST / CENTRAL
        my $GROUP3 = $split_info[6];    ## Human Feeding
        my $GROUP4 = $split_info[7];
        #my $TARGET = $split_info[0];

        $POPdir =~ s/\_CDS$//g;

        my ($country,$population) = split /\./, $POPdir;

        if ( $POPdir!~/aedes-albopictus.Jardin_Panteon_CDS/ && $POPdir!~/madagascar.Le_dauget_CDS/ ) {            
            # GROUP1 == FULL
            # ------------------------------------------------------------------------------------------------------------ 
            if ($GROUP1=~/Americas/ || $GROUP1=~/Asia/ || $GROUP1=~/Oceania/) {
                $MAJOR_POPULATION_GROUPS{Out_of_Africa}{$population} = $GROUP1; 
                #print ">>>$country\t$population\t|$POPdir|\n";
            }
            elsif ($GROUP1=~/Africa/) {
                #print "$line\n";
                $MAJOR_POPULATION_GROUPS{Africa}{$population} = $GROUP1;
            }
            else {
                #print "$line\n";
            }             
        }
    }
}








# read directories
# ---------------------------------------------------------------------------------------------------------------------------------------------------

my @main_directory_VCFs = qw(/path/to/results/selection_dnds/parsed_results/);

foreach my $i (@main_directory_VCFs){
    chdir $i;
    #print "[1] $i ... done\n";
    opendir DIR, ".";
    
    my @subdir= readdir DIR;
    &GoDown1(@subdir);
    closedir DIR;
}



# >>> GO-GO to NEXT SUBDIR >>>
sub GoDown1{
    my @subdir= ();
    @subdir= @_;

    if ($#subdir > 1){
        @subdir= @subdir[0..$#subdir];
        my @subdir2= sort {lc($a) cmp lc($b)} @subdir;

        my $POPdir = $main_DIR;
           $POPdir =~ s/\_CDS$//g;

        my ($country,$population) = split /\./, $POPdir;

        foreach my $j (@subdir2){
            chomp $j;
            
            # single consensus for a protein coding gene of a single population (~14,000 genes)
            # ------------------------------------
            if ( $j!~ /^\.$/ && $j!~ /^\.\.$/ && $j!~ /^\.listing$/ && $j=~/aedes_aegypti\.local_selection\.(.+)\.(.+)\.dnds_omega\.main_table\.missing_geneIDs\.txt$/ || $j=~/aedes_aegypti\.local_selection\.(.+)\.(.+).dnds_omega\.main_table\.txt$/ ) {
                #my ($targetfile,$country) = ("","");
                
                my $targetfile = $j;
                my $country    = $1;
                my $population = $2;

                my $MAJOR_GROUP = "";  my $MINOR_GROUP = "";
                
                if ( $MAJOR_POPULATION_GROUPS{Out_of_Africa}{$population} ) {  $MAJOR_GROUP = "Out_of_Africa";  }
                elsif ( $MAJOR_POPULATION_GROUPS{Africa}{$population} )     {  $MAJOR_GROUP = "Africa";  }
                else {  print "WARNING!\t$country\t$population\t$targetfile\n"  }
                #print "$MAJOR_GROUP\t$country\t$population\t$targetfile\n";


                open INFILE, "$targetfile" or die "CANNOT open INFILE $targetfile ~ $!\n";             
                

                # MAIN
                # ------------------------------------------------------------------------------------------
                while (<INFILE>) {
                    chomp;
                    my $line = $_;

                    my $selectionType_FLAG3;
                    my ($geneID,$Ka,$Ks,$KaKs,$Pval,$signif) = ("","","","","","");

                    if ($line!~/^\#.+/ && $line!~/^POP_TARGET.+/) {
                        my ($POP_TARGET,$COUNTRY,$POPULATION,$GENE_ID1,$GENE_ID2,$NON_SYNONYMOUS_N,$SYNONYMOUS_D,$OMEGA,$NON_SYNONYMOUS_N_RATIO,$SYNONYMOUS_D_RATIO,$KAPPA,$t_VALUE) = split /\t/, $line;

                        my $orthologous_Albopictus = "";
                        if ($HASH_ORTHOLOGS_GENES{$GENE_ID1}) {
                            $orthologous_Albopictus = $HASH_ORTHOLOGS_GENES{$GENE_ID1};
                        }
                        else {
                            print  "WARNING!!\t$MAJOR_GROUP\t$line\n";
                        }
                        
                        ## CUTOFF
                        # /// ------------------------------------------------------------------------------------------------------------------------------------ /// #
                        my $neutrality_T3_1A = 0.80;   my $neutrality_T3_1B = 1.20;

                        if ($OMEGA < 0.50 && $OMEGA!~/^NA$/)                                                 {  $selectionType_FLAG3 = "strong_negative";  }
                        elsif ($OMEGA >= 0.50 && $OMEGA < $neutrality_T3_1A && $OMEGA!~/^NA$/)               {  $selectionType_FLAG3 = "weak_negative";    }
                        elsif ($OMEGA >= $neutrality_T3_1A && $OMEGA <= $neutrality_T3_1B && $OMEGA!~/^NA$/) {  $selectionType_FLAG3 = "nearly_neutral";   }
                        elsif ($OMEGA >  $neutrality_T3_1B && $OMEGA <= 1.50 && $OMEGA!~/^NA$/)              {  $selectionType_FLAG3 = "weak_positive";    }
                        elsif ($OMEGA > 1.50 && $OMEGA!~/^NA$/)                                              {  $selectionType_FLAG3 = "strong_positive";  }
                        elsif ($OMEGA=~/^NA$/)                                                               {  $selectionType_FLAG3 = "ERROR"; }
                        else {  print "WARNING_THRESHOLD_3:\r$line\n";;  }
                        # /// ------------------------------------------------------------------------------------------------------------------------------------ /// #
                        



                        my $selection_type_SD7 = $selectionType_T3;
                           $selection_type_SD7 =~ s/\[//gis;
                           $selection_type_SD7 =~ s/\]//gis;
                           
                        print OUTFILE_SUPPLEMENTARY_DATA_SEVEN "$MAJOR_GROUP\t$COUNTRY\t$POPULATION\t$GENE_ID1\t$orthologous_Albopictus\$selectionType_FLAG3\t$NON_SYNONYMOUS_N\t$SYNONYMOUS_D\t$OMEGA\t$NON_SYNONYMOUS_N_RATIO\t$SYNONYMOUS_D_RATIO\t$KAPPA\t$t_VALUE\n";

                }
                close INFILE;
                # ------------------------------------------------------------------------------------------
            }
        }
    }
}
close OUTFILE_GENE_DNDS_SELECTION;
close OUTFILE_SUPPLEMENTARY_DATA_SEVEN;



