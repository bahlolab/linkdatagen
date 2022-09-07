# linkdatagen

Linkdatagen website: https://bioinf.wehi.edu.au/software/linkdatagen/index.html

LINKDATAGEN is a PERL (5.8+) script that generates datasets for linkage analysis, relatedness checking, IBD and HBD inference and association analyses using genotypes generated from any of Affymetrix SNP chips, Illumina SNP chips or massively parallel sequencing (MPS; otherwise known as next generation sequencing or NGS) data.

LINKDATAGEN selects markers for linkage and association analyses and performs nuclear family Mendelian error detection. It also performs sex checks using both X and Y chromosome markers (if available). Linkdatagen supports all eleven HapMap populations (ASW, CEU, CHB, CHD, GIH, JPT, LWK, MEX, MKK, TSI and YRI) for the following platforms, unless otherwise specified:

-   Affymetrix: 50K Xba, 50K Hind, 250K Sty, 250K Nsp, 5.0, 6.0 and 500K Sty+Nsp
-   Illumina: 610Quad, Cyto12, Omni Express, 1M, plus others due to high overlap between Illumina platforms
-   MPS data: Hapmap2 (for CEU, CHB, JPT and YRI populations) and HapMap3
-   PLINK data (TODO: should this be included here?)

However, many other chips can be catered for (see -chip section in manual). 

LINKDATAGEN creates output files for the linkage mapping software: ALLEGRO, MERLIN, MORGAN and PLINK, as well as for BEAGLE, FESTIM, PREST, fastPHASE and RELATE. Many of these programs are available through the Rockefeller website (<http://linkage.rockefeller.edu/soft>).

# Cite

If you use LINKDATAGEN, please acknowledge by citing:

> Bahlo M, Bromhead CJ. Generating linkage mapping files from Affymetrix SNP chip data. *Bioinformatics* 2009;25(15):1961-2.

If you use the LINKDATAGEN -data m option along with vcf2linkdatagen.pl then please also cite:

> Smith KR, Bromhead CJ, Hildebrand MS, Shearer AE, Lockhart PJ, Najmabadi H, Leventer RJ, McGillivray G, Amor DJ, Smith RJ, Bahlo M. Reducing the exome search space for Mendelian diseases using genetic linkage analysis of exome genotypes. *Genome Biology* 2011;12:R85.


# Developer notes

## git branch structure

- main: Scripts ready for public use.
- beta: Scripts for testing within the Bahlo lab or when you have been instructed to use the beta version.
- (others): Development branches, not ready for production.
