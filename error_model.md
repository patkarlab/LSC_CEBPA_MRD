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