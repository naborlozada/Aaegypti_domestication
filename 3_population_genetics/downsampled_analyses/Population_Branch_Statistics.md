
## Population Branch Statistics (PBS)

After calculating the SAF for each population, as explained in this [thread](https://github.com/naborlozada/Aaegypti_domestication/blob/main/downsampled_analyses/Genetic_diversity_metrics.md), we calculated the SFS, the FST, and the PBS based on populations pairwise comparisons as follows.

First, we created groups of populations that represent three focal populations:

A) Target focal branch population: Out of Africa (OoA)
B) Outgroup focal branch population: Eastern Africa (EAS) without human feeding domesticated Aaa-like mosquitoes from Rabai, Kenya (RABd).
C) Sister focal branch population: Western Africa (WES), which was divided in two subgroups of datasets: 
     1) Western Africa without the human feeding mosquitoes THI and NGY from Senegal (WESxHF), and 
     2) Western Africa with all samples, including the human feeding Aaa-like mosquitoes (WESF).
D) Sister focal branch population: The Aaa-like human feeding mosquitoes (HF): RABd, THI and NGY.

Comparative phylogenetic distance analyses:
* (B;C1,A)
* (B;C2,A)
* (B;D,A)



```bash

# case: (B;C1,A)
# ----------------------------------------------------------

# Calculate the 2Dsfs prior
realSFS popA.saf.idx popB.saf.idx > popA_popB.sfs;
realSFS popA.saf.idx popC1.saf.idx > popA_popC1.sfs;
realSFS popB.saf.idx popC1.saf.idx > popB_popC1.sfs;


# ----------------------------------------------------------
# Case: B; A vs C1
# ----------------------------------------------------------

# Prepare the fst for persite/window analysis
realSFS fst index popA.saf.idx popB.saf.idx popC1.saf.idx -sfs popA_popB.sfs -sfs popA_popC1.sfs -sfs popB_popC1.sfs -fstout popA_popB_popC1.fst

# Get the global estimate
realSFS fst stats popA_popB_popC1.fst.idx > popA_popB_popC1.pbs.global_stats.txt 2>&1 | tee popA_popB_popC1.pbs.global_stats.stderr.log

# Sliding windows 100kb
realSFS fst stats2 popA_popB_popC1.fst.idx -win 99999 -step 100000 -type 2 > popA_popB_popC1.sliwin_100kb.pbs 2>&1 | tee popA_popB_popC1.pbs.sliwin_100kb.stderr.log



# ----------------------------------------------------------
# Case: B; A vs C2
# ----------------------------------------------------------

# Prepare the fst for persite/window analysis
realSFS fst index popA.saf.idx popB.saf.idx popC2.saf.idx -sfs popA_popB.sfs -sfs popA_popC2.sfs -sfs popB_popC2.sfs -fstout popA_popB_popC2.fst

# Get the global estimate
realSFS fst stats popA_popB_popC2.fst.idx > popA_popB_popC2.pbs.global_stats.txt 2>&1 | tee popA_popB_popC2.pbs.global_stats.stderr.log

# Sliding windows 100kb
realSFS fst stats2 popA_popB_popC2.fst.idx -win 99999 -step 100000 -type 2 > popA_popB_popC2.sliwin_100kb.pbs 2>&1 | tee popA_popB_popC2.pbs.sliwin_100kb.stderr.log




# ----------------------------------------------------------
# Case: B; A vs D
# ----------------------------------------------------------

# Prepare the fst for persite/window analysis
realSFS fst index popA.saf.idx popB.saf.idx popD.saf.idx -sfs popA_popB.sfs -sfs popA_popD.sfs -sfs popB_popD.sfs -fstout popA_popB_popD.fst

# Get the global estimate
realSFS fst stats popA_popB_popD.fst.idx > popA_popB_popD.pbs.global_stats.txt 2>&1 | tee popA_popB_popD.pbs.global_stats.stderr.log

# Sliding windows 100kb
realSFS fst stats2 popA_popB_popD.fst.idx -win 99999 -step 100000 -type 2 > popA_popB_popD.sliwin_100kb.pbs 2>&1 | tee popA_popB_popD.pbs.sliwin_100kb.stderr.log


```
