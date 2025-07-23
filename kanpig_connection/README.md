## Kanpig connector from Liza's scripts

Main script: `listing.sh`. 
Dependencies:
```
bedgovcf
htslib
samtools
bcftools
bedtools
kanpig
```

Takes 3 arguments:
1. `.bed` file location from fly catcher;
2. `.bam` file location for the aforementioned `.bed` file;
3. `.yaml` config for `bedgovcf` (provided as `bed2vcf.yaml` here).

Results will be in `fc_kanpig_results`. Temporary directory is also provided w/ `.bam` indexes created by kanpig for future use.

Workflow:
1. Fetches sequences from reference genome provided in `support` (now GRCh38, subject to change);
2. Converts `.bed` to `.vcf` with some tags (GT field is imputed as ./1 rn. Once we have a way to estimate dosage, this can change);
3. Runs `kanpig` and filters entries with empty GT field;

Needs:
1. GRCh38 indexed in `support`;
2. ??? Nothing else I guess.
