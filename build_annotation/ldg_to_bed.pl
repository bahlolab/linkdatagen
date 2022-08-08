# Convert from Linkdatagen universal format to a BED file with columns collapsedin the BED 'name' field.
#

use 5.020;
use strict;
use warnings;


while(<>) {
    if(/^#/ or /^Chrom/) {
        next;
    }
    my @line = split /\t/;
    say $line[0] . "\t" . ($line[1] - 1) . "\t$line[1]\t$line[2]\t1\t$line[3]";
}

