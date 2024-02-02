
## RAiSD detection of hard selective sweeps.

The detection of the selective sweeps was performed at two levels: local (population) and global level.

## **LOCAL**
For each of the 40 populations, RAiSD was applied with the following commandline structure:

```bash
myPOPS=/home/nlozada/selection/vcf_byPopulation/*.vcf.gz

for i in $myPOPS; do
    POP=$(echo $i | cut -d'.' -f 1);
    if [[ $i =~ chrm1 ]]; then
        RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm1.raisd.win50 -I $i;
        wait; sleep 2;
    elif [[ $i =~ chrm2 ]]; then
        RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm2.raisd.win50 -I $i;
        wait; sleep 2;
    elif [[ $i =~ chrm3 ]]; then
        RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm3.raisd.win50 -I $i; 
    else
        echo "WARNING FILE!" $i;
    fi
done

# Create for each chrosomosome a RAiSD outout, settings:
#    * [-y 2] doploidy 
#    * [-A 0.905] do an internal test and set a probability value to be used for the quantile function in R, and generates a Manhattan plot for the final mu-statistic score using an internal Rscript.
#    * [-M 0] work with NO missing data
#    * [-w 50] use sliding window of 50
#    * [-s, -R, -P, -D, -O] generate a report file, include additional information (window start and end, and the mu-statistic factors for variation, SFS, and LD), genetarate internal exploratory plots, show job progress. 
```

## **GLOBAL**

All populations were grouped in three major groups: 
1) Out of Africa
2) Africa (all populations)
3) Africa (all populations, except for human feeding mosquitoes RABd, THI, NGY)

For each group, the same RAiSD commandline from above was applied separately.

```bash

# vcf merged all pops
mainVCF=aaegL5.all_populations.merged.biallelic.geno80.vcf.gz;

# Out of Africa:
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm1.raisd.win50 -I mainVCF -S popList.out_of_africa.txt;
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm2.raisd.win50 -I mainVCF -S popList.out_of_africa.txt;
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm3.raisd.win50 -I mainVCF -S popList.out_of_africa.txt;
# africa
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm1.raisd.win50 -I mainVCF -S popList.africa_full.txt;
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm2.raisd.win50 -I mainVCF -S popList.africa_full.txt;
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm3.raisd.win50 -I mainVCF -S popList.africa_full.txt;
# human feeding
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm1.raisd.win50 -I mainVCF -S popList.human_feeding.txt;
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm2.raisd.win50 -I mainVCF -S popList.human_feeding.txt;
RAiSD -s -R -P -D -O -A 0.905 -M 0 -y 2 -w 50 -n ${POP}.geno80.chrm3.raisd.win50 -I mainVCF -S popList.human_feeding.txt;
```

use the R script `outliers_RAiSD.R` to parse and obtain outliers with statistical support.
