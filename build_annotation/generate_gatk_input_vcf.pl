# Filter VCF file to just the SNPs we have annotated.

use 5.020;
use strict;
use warnings;
use autodie;

open my $annot, '<', 'annotHapMap3U_hg38.txt';
my %snps;

while(<$annot>) {
    my @line = split /\t/;
    $snps{$line[2]} = 1;
}
close $annot;

open my $dbsnp, '<', 'data/dbsnp_138.hg38.vcf';

while(<$dbsnp>) {
    if(/^#/) {
        print;
        next;
    }
    my @line = split /\t/;
    if(defined($snps{$line[2]})) {
        print;
    }
}
