# Generation of Error model using mean background error
Calculation of mean background error in this study was done as described in Thol F et. al. Blood (2018) 132 (16): 1703–1713.
[PMID: 30190321](https://doi.org/10.1182/blood-2018-02-829911 "doi link")

## Steps for generating the background error rate

1. Median calculation across callers

    - For each variant, median REF count, ALT count and VAF% are calculated across the three callers.

3. LVAF determination (Largest Variant Allele Frequency)

    - Variants are sorted by chromosome and position, and the highest median VAF per genomic position is retained as LVAF%.
        
    - Only SNPs and short indels (≤2 bp length difference) are used for background modeling.

4. Background error estimation

    LVAF% values are used to calculate:

    - Mean background error
    - Standard deviation

    Outliers beyond mean ± 2.5 SD are removed before final estimation.

5. MRD Threshold

    MRD threshold is defined as:

    - Background + 3 × SD
    - Variants exceeding this threshold are labelled PASS; others are labelled FAIL.
    - Large indels are labelled NA.

Thresholds are defined as described in the supplementary methods of the [article](https://doi.org/10.1182/blood-2018-02-829911 "doi link").
