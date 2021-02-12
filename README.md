# AneuploidyScore
The point of this package is to calculate the **Aneuploidy Score** for a chromosomal arm SCNA/aneuploidy (**CAA**) as detailed in the following two papers:

  - Shukla, A., Nguyen, T. H., Moka, S. B., Ellis, J. J., Grady, J. P., Oey, H., ... & Duijf, P. H. (2020). Chromosome arm aneuploidies shape tumour evolution and drug response. Nature communications, 11(1), 1-14. (https://doi.org/10.1038/s41467-020-14286-0)
 
  - Cohen-Sharir, Y., McFarland, J. M., Abdusamad, M., Marquis, C., Bernhard, S. V., Kazachkova, M., ... & Ben-David, U. (2021). Aneuploidy renders cancer cells vulnerable to mitotic checkpoint inhibition. Nature, 1-6. (https://doi.org/10.1038/s41586-020-03114-6)
  
R package to calculate the Aneuploidy Score from Chromosome Arm-level SCNAs/Aneuploidies (CAAs) as outlined and expanded by Shukla et al. (https://doi.org/10.1038/s41467-020-14286-0)

# Installation
Install using
```
devtools::install_github("quevedor2/aneuploidy_score", ref='master') # Latest stable build
devtools::install_github("quevedor2/aneuploidy_score", ref='dev') # Development branch
```

A guide for how to use the package can be found at `inst/rmd/demo_guide.html`

# Details of Aneuploidy Score/CAA  
Defintions:
 
  - Aneuploidy score (**AS**): as described by Cohen-Sharir et al., is simply the total number of arm-level gains and losses for a tumor, adjusted for ploidy. 
 
  - Chromosomal arm-level SCNAs/Aneuploidies (**CAA**): as described by Shukla et al., are arm-level aneuploidies that are pervasive in cancer and found ~30x higher than expected based on the inverse-length distribution of focal SCNAs.


AS and CAAs are calculated two different ways:

  - Shukla method: 
    - A **Log 2 Ratio (L2R)**-based method
    - Identify gains and losses using a 0.2 and -0.2 threshold respectively
    - If the CN segment intersects the centromere positions, discard it (Line 118-126 (commit: 25ea7bd): https://github.com/pascalduijf/CAAs_1/blob/master/scripts/CAAs_from_cell_lines.py)
    - Flag the arm as a CAA if the total length for the gain/loss segment in the chromosome arm exceeds **90% (0.90)** of the total chromosome arm
    - [Optional]: whole-genome doubling is determined if half of the autosomal genome had two or more copies of the more frequent (maternal or paternal) allele.
  
  - Cohen-Sharir method: 
    - A **Total Copy Number (TCN)**-based method
    - If a CN-segment spans the centromere, split it and assign the segment to its respective arm.
    - Take the segment-size weighted median total copy number of each chromosome arm (*arm_wMedian*)
    - Assign the *arm_wMedian* a loss/neutral/gain (-1/0/1) status based on relation to the cell lines *ploidy*
    - Using a 0.5 threshold (from the default round() function), determine whether the *arm_wMedian* is >, =, or < than the *ploidy*
    - **AS** = Total number of arm gains/losses
    


