# NCA_2nd Round

## Data by age
0. Total number of samples: 20 samples --> 3 invalid --> **17 valid samples**
1. 3 weeks: 6 samples
2. 6 weeks: 6 samples --> 2 of them are invalid --> 4 valid samples
3. 9 weeks: 4 samples --> 1 of them is invalid --> 3 valid samples
4. 12 weeks: 4 samples




|Age | Sample line    | Date     | Tech  | Index           | df name  |Space|Remarks
|---:|:--------------:|:--------:|:-----:|:---------------:|:--------:|:---:|:------:|
|6W  |45F_01          |8/4/2020  | indrop| 1 - CTCTCTAT    |w6_45f.01 |whole|invalide|
|6W  |47F_01          |8/4/2020  | indrop| 2 - TATCCTCT    |w6_47f.01 |whole|invalide|
|6W  |G10             |8/4/2020  | indrop| 3 - GTAAGGAG    |w6_g10    |whole|        |
|6W  |WT              |8/4/2020  | indrop| 4 - ACTGCATA    |w6_wt     |whole|        |
|6W  |45F_02          |8/4/2020  | indrop| 5 - AAGGAGTA    |w6_45f.02 |whole|        |
|6W  |47F_02          |8/4/2020  | indrop| 6 - CTAAGCCT    |w6_47f.02 |whole|        |
|3W  |45F             |8/21/2020 | indrop| 7 - CGTCTAAT    |w3_45f    |whole|        |
|3W  |47F             |8/21/2020 | indrop| 8 - TCTCTCCG    |w3_47f    |whole|        |
|3W  |G10             |8/21/2020 | indrop| 9 - TCGACTAG    |w3_g10    |whole|        |
|3W  |WT              |8/21/2020 | indrop| 10 - TTCTAGCT   |w3_wt     |whole|        |
|3W  |45F nuc         |8/21/2020 | indrop| 11 - CCTAGAGT   |w3_45f_nuc|nuc  |        |
|3W  |47F nuc         |8/21/2020 | indrop| 12 - GCGTAAGA   |w3_47f_nuc|nuc  |        |
|12w |45f             |10/2/2020 | indrop| 13 - CTTAATAG   |w12_45f   |whole|        |
|12w |47f             |10/2/2020 | indrop| 14 - ATAGCCTT   |w12_47f   |whole|        |
|12w |g10             |10/2/2020 | indrop| 15 - TAAGGCTC   |w12_g10   |whole|        |
|12w |wt              |10/2/2020 | indrop| 16 - TCGCATAA   |w12_wt    |whole|        |
|9w  |45f             |10/8/2020 | indrop| 17 - TTACCTCC   |w9_45f    |whole|invalide|
|9w  |47f             |10/8/2020 | indrop| 18 - CAGTTATG   |w9_47f    |whole|        |
|9w  |g10             |10/8/2020 | indrop| 19 - CCTTTACT   |w9_g10    |whole|        |
|9w  |wt              |10/8/2020 | indrop| 21 - GAGACGGA   |w9_wt     |whole|        |


Invalid data files: w6_45f.01 and w6_47f.01
Samples included in the downstream analysis: 18 samples 

* 3w 4 samples
* 3w 2 nuc seqs (45f & 47f)
* 6w 4 samples
* 9w 4 samples
* 12w 4 samples


## Quality metrics

**1. Cell counts by sample**

<img src="figures/QC/QC metrics/ncells.jpeg" width="800">
 
**2. Number of UMIs/transcripts per cell**

<img src="figures/QC/QC metrics/umi_counts.jpeg" width="800">

**3. Detected unique genes per cell**
<img src="figures/QC/QC metrics/unique_genes_box.jpeg" width = "1000">

<img src="figures/QC/QC metrics/unique_genes_box_log.jpeg" width = "1000">

**4. UMI-Unique genes-Mitochondrial content**
<img src="figures/QC/QC metrics/umi_genes_mt.jpeg" width = "1200">

**5. Complexity**
<img src="figures/QC/QC metrics/complexity.jpeg" width = "800">
