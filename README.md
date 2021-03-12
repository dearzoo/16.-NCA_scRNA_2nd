# NCA_2nd Round

## Data by age
##### Start with **17 valid samples**
1. 3 weeks: 6 samples
2. 6 weeks: 6 samples --> 2 of them are invalid --> 4 valid samples
3. 9 weeks: 4 samples --> 1 of them is invalid --> 3 valid samples
4. 12 weeks: 4 samples



| No. | Age | Sample line |    Date   |  Tech  |     Index     |   df name  | Space |       Remarks       |
|:---:|:---:|:-----------:|:---------:|:------:|:-------------:|:----------:|:-----:|:-------------------:|
|  1  |  6W |    45F_01   |  8/4/2020 | indrop |  1 - CTCTCTAT |  w6_45f.01 | whole |     **invalid data file**     |
|  2  |  6W |    47F_01   |  8/4/2020 | indrop |  2 - TATCCTCT |  w6_47f.01 | whole |     **invalid data file**     |
|  3  |  6W |     G10     |  8/4/2020 | indrop |  3 - GTAAGGAG |   w6_g10   | whole |                     |
|  4  |  6W |      WT     |  8/4/2020 | indrop |  4 - ACTGCATA |    w6_wt   | whole |                     |
|  5  |  6W |    45F_02   |  8/4/2020 | indrop |  5 - AAGGAGTA |  w6_45f.02 | whole |                     |
|  6  |  6W |    47F_02   |  8/4/2020 | indrop |  6 - CTAAGCCT |  w6_47f.02 | whole |                     |
|  7  |  3W |     45F     | 8/21/2020 | indrop |  7 - CGTCTAAT |   w3_45f   | whole |                     |
|  8  |  3W |     47F     | 8/21/2020 | indrop |  8 - TCTCTCCG |   w3_47f   | whole |                     |
|  9  |  3W |     G10     | 8/21/2020 | indrop |  9 - TCGACTAG |   w3_g10   | whole |                     |
|  10 |  3W |      WT     | 8/21/2020 | indrop | 10 - TTCTAGCT |    w3_wt   | whole |                     |
|  11 |  3W |   45F nuc   | 8/21/2020 | indrop | 11 - CCTAGAGT | w3_45f_nuc |  nuc  |                     |
|  12 |  3W |   47F nuc   | 8/21/2020 | indrop | 12 - GCGTAAGA | w3_47f_nuc |  nuc  |                     |
|  13 | 12w |     45f     | 10/2/2020 | indrop | 13 - CTTAATAG |   w12_45f  | whole |                     |
|  14 | 12w |     47f     | 10/2/2020 | indrop | 14 - ATAGCCTT |   w12_47f  | whole |                     |
|  15 | 12w |     g10     | 10/2/2020 | indrop | 15 - TAAGGCTC |   w12_g10  | whole |                     |
|  16 | 12w |      wt     | 10/2/2020 | indrop | 16 - TCGCATAA |   w12_wt   | whole |                     |
|  17 |  9w |     45f     | 10/8/2020 | indrop | 17 - TTACCTCC |   w9_45f   | whole | **gene counts too low** |
|  18 |  9w |     47f     | 10/8/2020 | indrop | 18 - CAGTTATG |   w9_47f   | whole |                     |
|  19 |  9w |     g10     | 10/8/2020 | indrop | 19 - CCTTTACT |   w9_g10   | whole |                     |
|  20 |  9w |      wt     | 10/8/2020 | indrop | 21 - GAGACGGA |    w9_wt   | whole |                     |



## Quality metrics (17 datasets)

**1. Cell counts by sample**

<img src="figures/QC/QC metrics/ncells.jpeg" width="800">
 
**2. Number of UMIs/transcripts per cell**

<img src="figures/QC/QC metrics/umi_counts.jpeg" width="800">

**3. Detected unique genes per cell**

- The 12-week g10 sample shows a way lower gene counts than other samples

<img src="figures/QC/QC metrics/unique_genes_box.jpeg" width = "1000">

<img src="figures/QC/QC metrics/unique_genes_box_log.jpeg" width = "1000">


**Conclusion**
 - better to exclude w12_g10

 


- **Median value of unique gene counts by sample**

|   Sample   | Median_unique gene counts |
|:----------:|:-------------------------:|
|   w6_g10   |           326.0           |
|    w6_wt   |           301.0           |
|  w6_45f.02 |           402.0           |
|  w6_47f.02 |           388.0           |
|   w3_45f   |           258.0           |
|   w3_47f   |           337.0           |
|   w3_g10   |           275.0           |
|    w3_wt   |           313.0           |
| w3_45f_nuc |           238.0           |
| w3_47f_nuc |           235.0           |
|   w12_45f  |           264.5           |
|   w12_47f  |           244.5           |
|   w12_g10  |            11.5           |
|   w12_wt   |           223.0           |
|   w9_45f   |            0.0            |
|   w9_47f   |           128.0           |
|   w9_g10   |           219.0           |
|    w9_wt   |           242.5           |

**4. UMI-Unique genes-Mitochondrial content**

- The 12-week g10 sample shows low quality characteristics. 
- Positioned under the set intercepts.

<img src="figures/QC/QC metrics/umi_genes_mt.jpeg" width = "1200">

**5. Complexity**

- The 12-week g10 sample shows very low complexity when compared with other samples. 
<img src="figures/QC/QC metrics/complexity.jpeg" width = "800">


## Violin plots by sample

- The 3rd columns are mitrochondrial content %. (sample names were added only for the convenience purpose)
- [Violin Plots](https://github.com/dearzoo/16.-NCA_scRNA_2nd/blob/main/figures/violin/qc_violin_with_names.pdf) for the violin plots



#### Sample-level filtration
**12-week 45 sample** looks low quality --> will be additionally ***excluded*** from downstream analysis





