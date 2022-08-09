#!/usr/bin/perl -w

my $script_revision = "3.1.0";
my $last_changed_date = "2022-08-09 16:50 AEST";

my $license = <<LICENSE;
    Copyright (C) 2009-2022 Melanie Bahlo, Catherine Bromhead, Thomas Scerri,
    Katherine Smith, Rick Tankard and Luke Gandolfo. 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
LICENSE

#NOTE: This is the *new* linkdatagen which amalgamates linkdatagen_affy, linkdatagen_illumina and linkdatagen_mps
# not to be confused with the original linkdatagen which only dealt with Affy data

=head1 LINKDATAGEN

                    ------------------------------- LINKDATAGEN -------------------------------
                    now incorporating:
                    LINKDATAGEN_AFFY
                    LINKDATAGEN_ILLUMINA
                    LINKDATAGEN_MPS
                    M Bahlo, CJ Bromhead
        
                    with help from Tom Scerri, Katherine Smith, Rick Tankard and Luke Gandolfo.
                    ---------------------------------------------------------------------------

    If you use LINKDATAGEN please acknowledge by citing:
    Bahlo M, Bromhead CJ. Generating linkage mapping files from Affymetrix SNP chip data. Bioinformatics 2009;25(15):1961-2.
    
    If you use LINKDATAGEN -data m (for MPS data) please acknowledge by citing:
    Smith KR, Bromhead CJ, Hildebrand MS, Shearer AE, Lockhart PJ, Najmabadi H, Leventer RJ, McGillivray G, Amor DJ, Smith RJ, Bahlo M. 
    Reducing the exome search space for Mendelian diseases using genetic linkage analysis of exome genotypes. Genome Biology 2011;12:R85.

=head1 License

This script is free software that is licensed under GPL v2 and comes with no warranty. For details, run the -license option. 

=head1 Help

Please see documentation available from http://bioinf.wehi.edu/software/linkdatagen/ for more detailed help & instructions.

All input files and output files may be gzip or bzip2 compressed when the file name is given (with extension gz or bz2), provided the respective module PerlIO::gzip or 
PerlIO::via::Bzip2 is installed from CPAN. 

  <---------- ---------- --------- ---------- ---------- ---------- ----------    Increase window size until this is on one line    ---------- ---------- ---------- ---------- ---------- --------->

=head2 Mandatory Options:

  -annotFile < filename >
  -chip { 1 , 2 , 3 , 4 , 5 , 6 , 7 }

      Both -annotFile and -chip are used to specify an annotation file.  They are mutually exclusive; only one must be specified.  If -chip is used, then -annotDir must also be used.
      -annotFile directly specifies an annotation file.
      -chip dictates the annotation file required relative to -data.

      Affymetrix SNP chip data has 7 possible values: 1 = 50k_Xba, 2 = 50k_Hind, 3 = 250k_Sty, 4 = 250k_Nsp, 5 = 5.0, 6 = 6.0, 7 = 500k_Nsp+Sty
      Illumina SNP chip data has 6 possible values: 1 = 370Duo, 2 = 610Quad, 3 = 660Quad, 4 = Cyto12, 5 = Omni_Express, 6 = 1M
      MPS data has 2 possible values: 1 = HapMap2 (only CEU, YRI, HAN, JPT, ~4 million SNPs), 2 = HapMap3 (all eleven HAPMAP pops, ~1.5 million SNPs)
  
  -callDir < path to directory >
  -callFile < filename >

      Both -callDir and -callFile are used to locate the genotype data.  They are mutually exclusive; only one must be specified.
      -callDir should be used with Illumina SNP chip data and specifies the location of the genotyping data files, named "*_Final_Report.txt".
      -callFile specifies a single file containing genotype data in either BRLMM (typically brlmm.calls.txt) or CRLMM (typically crlmm.calls.txt) formats.
      -callFile is to be used with Affymetrix SNP chip data (both BRLMM and CRLMM) and MPS data (BRLMM only).
  
  -data { a, i, m }

      This indicates the source of the genotype data, where a = Affymetrix, i = Illumina, m = MPS (or NGS). Example: -data a for Affymetrix SNP chip data.
  
  -pedfile < filename >

      The pedigree file with 6 columns: pedID/indID/father's_ID/mother's_ID/sex/affection_status.
  
  -freq
  -prog { all | al | me | mo | pl | pr | cp | be | fe | re | fp }
  -popHetTest { summary , verbose , perChr , perChrVerbose }
  -bestPopTest

      At least one of EITHER -freq and/or -prog and/or -popHetTest must be selected, OR otherwise -bestPopTest alone. 
  
    -freq will output the allele frequency estimates from the founders of the given samples' genotype calls.
      -prog specifies the program choice(s) for output file formats in a comma-seperated list, e.g. -prog me,pl,re,mo, where:
  
        fp = fastPHASE, pr = PREST, al = ALLEGRO, me = MERLIN, mo = MORGAN, pl = PLINK, be = BEAGLE, fe = FESTIM, re = RELATE, cp = internal format, all = all formats
  
    -popHetTest runs a goodness-of-fit test by comparing the expected allele frequencies of the selected population against the observed allele frequencies in the genotype call files.
      -bestPopTest runs the same goodness-of-fit as -popHetTest using all available population allele frequencies, and reports the best fitting population for each sample.
  
  -whichSamplesFile < filename >
  -whichSamplesList < filename >

      Both -whichSamplesFile and -whichSamplesFile are used to link the individuals' genotype data with their position in the pedfile.  They are mutually exclusive; only one must be specified.
      -whichSamplesFile specifies a file containing a single line of N space-delimited values, where N equals the number of individuals in the pedfile.
      Each value corresponds to an individual in the pedfile (in the same order as the pedfile), where 0 indicates no genotype data available for that individual.
      Otherwise the value is one of {1 , ... , genotyped samples} where each value indicates a column in the BRLMM file.
      -whichSamplesList specifies a file containing one line for each individual in the pedfile, in the same order as the pedfile.
      Each line should specify the path to the genotyping file for that individual.  The line should be 0 if there is no sample available.
      If there are pedigree errors then it is worthwhile double-checking either the -whichSamplesFile or -whichSamplesList file.
  
  -MendelErrors { keep | removeSNPs | removeGPF }

      Specify action to be taken with SNP markers displaying Mendelian errors.
      "keep" - all SNPs (and associated genotypes) are kept regardless of Mendelian errors.
      "removeSNPs" - SNPs (and all associated genotypes) displaying a Mendelian error in any familiy are removed from the entire dataset (of all families).
      "removeGPF" - "genotypes per family" are wiped ("made missing") for any SNPs displaying Mendelian errors.  This occurs only for that specific family.
  
=head2 Other Options:
  
  -actg

      Not functioning yet.
  
  -annotDir < path to directory >

      Specifies the location of the annotation file(s).  If using -chip (see above), then it must be used.
  
  -binsize { real number >= 0.0 }

      Specify size of bins, in cM, from which to choose each SNP. Default = 0.3 cM. Set "-binsize 0" will include all SNPs.
  
  -crlmm

      Specify data is in CRLMM format, where genotypes are coded 0 (missing), 1, 2, 3.  BRLMM data is coded -1 (missing), 0, 1 or 3.
      -crlmm is only relevant for "-data a" and is unlikely to be needed.
  
  -fileKeepSNPs < filename >

      Specify a list of SNPs by name (e.g. rs123456), one-per-line, that LINKDATAGEN may subsequently select from; i.e. it is an inclusion list.
  
  -fileRemoveSNPs < filename >

      Specify a list of SNPs by name (e.g. rs123456), one-per-line, that LINKDATAGEN may NOT subsequently select from; i.e. it is an exclusion list.
      SNPs listed in a -fileRemoveSNPs file will be excluded even if they are in the -fileKeepSNPs file, without warning.
  
  -help

      Print this help page to screen.
  
  -minDist { real number >= 0.0 }

      Specify the minimum distance (cM) allowable between selected SNPs.
      Default behaviour is SNPs must at least 0.2 cM, or 0.5 * binsize (cM), whichever is the smaller.
  
  -noX

      Declare that all chromosome X data to be excluded from analysis and/or output.
  
  -outputDir < prefix >

      Specify a prefix to the name of directories into which all files will be written.
      A directory called "<prefix>_tables" and containing log files will be created or written.
      Other directories called "<prefix>_<prog>" may also be created or written and will contain the files relevant to each program.
  
  -pop { ASW , CEU , CHB , CHD , GIH , JPT , LWK , MEX , MKK , TSI , YRI }
  -popCol { integer >= 1 }

      Both options specify a column of population allele frequenices from an annotation file.  They are mutually exclusive; only one may be specified.
      -pop is the most likely option to use as it is refers to the annotation files create by us.
      The default behaviour of LINKDATAGEN is to use the CEU population allele frequenices.  The different possible values for -pop correspond to:
  
          ASW : African ancestry in Southwest USA
          CEU : Utah residents with Northern and Western European ancestry from the CEPH collection
          CHB : Han Chinese in Beijing, China
          CHD : Chinese in Metropolitan Denver, Colorado
          GIH : Gujarati Indians in Houston, Texas
          JPT : Japanese in Tokyo, Japan
          LWK : Luhya in Webuye, Kenya
          MEX : Mexican ancestry in Los Angeles, California
          MKK : Maasai in Kinyawa, Kenya
          TSI : Toscans in Italy
          YRI : Yoruba in Ibadan, Nigeria (West Africa)
  
    -popCol should be used to select column K from a custom-made annotation file that is similar in structure to our annotation files, where the Kth column contains the desired allele frequencies.
  
  -regions { #,chr#,#:####-####,chr#:####-####,... }

      Specify comma-separated whole chromosomes and/or regions of chromosomes for subsequent analysis.
      Format for a whole chromosome is either "#" or "chr#" (that is, the "chr" is optional).  Here, the # must be a number from 1 to 22, or else the letter X, Y or M.
      Format for a region is "chr#:#######-#######" ("chr" is optional again).  Here, #######-####### are integers representing start and stop positions (bp).
      For example, "-regions chr5,X:12000-45000,3,chrY:3400-4300,chrM:100-200,6,chr19:7000-10000,chr19:30000-500000".
  
  -regionsFile < filename >

      Specify a file containing whole chromosomes and/or regions of chromosomes, one-per-line, for subsequent analysis.
      Format is at least identical to the -regions option (see above), although the use of ":" and "," to specify regions is optional.  For example, a file might contain:
  
        chr5
          3
          chrX:12000-45000
          chrY 3400 4300
          chrM:100-200
  
  -removeAIS

      Remove SNP markers that are absent from any Illumina input file. To be used with caution. Please firtsly understand why some markers are absent before continuing to use this option.
  
  -removeWFHBS { i , u }

      Remove "within-family homozygosity-by-state" SNP markers. To be used with caution. Please refer to the manual regarding use.
      When handling multiple families, the values "i" (for intersection) and "u" (for union) will remove the intersection or union of WFHBS markers across families, respectively.
      The default when applying -removeWFHBS is to remove the union, that is any SNP marker showing "within-family homozygosity-by-state" in any family.
  
  -randomSNP

      Declare that a random SNP should be selected from each bin interval.
      Default behaviour is to select the most heterozygous SNP based on the selected population allele frequencies.
  
=head2 Basic example:
  
  linkdatagen.pl -data a -pedfile ped.txt -whichSamplesFile wsf.txt -callFile genotypes.txt -chip 3 -annotDir ../annotationFiles/ -prog me,pl,fe
  
=cut 

use 5.006;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

# Compressed file opening modules, allows script to work even if modules are not installed:
my $perliogzip = 0;
eval { require PerlIO::gzip };
unless($@) { $perliogzip = 1 };
my $perlioviabzip2 = 0;
eval { require PerlIO::via::Bzip2 };
unless($@) { $perlioviabzip2 = 1 };
# predeclare subs
sub copen;
sub print_usage;

my $data;
my $randomSNP;
my $pedfile;
my $annotFile;
my $annotFileIlluminaType = 'univ';
my $whichSamplesFile;
my $whichSamplesList; #CATHEDIT
my $callDir;
my $callFile;
my $chip;
my $binsize;
my $outputDir;
my $annotDir;
my $crlmm;
my $brlmm;
my @prog;
my %progs = ("al"=>0,"me"=>0,"pr"=>0,"cp"=>0,"pl"=>0,"mo"=>0,"be"=>0,"fe"=>0,"re"=>0,"fp"=>0,"all"=>0);
#my $keepME;
my $MEO;	#Mendelian Error Option
my $removeWFHBS;
my $removeAIS;
my $minDist;
my $regions;
my $regionsFile;
my $pop;
my $popCol;
my $freq;
my $tmp;
my $noX;  #CATHEDIT
my $numchrom;  #CATHEDIT
my $start_time=0;
my $end_time=0;
my $popHetTest;
my $actg;
my $fileKeepSNPs;
my $fileRemoveSNPs;
my $seed;
my $help;
my $bestPopTest;
my $plSim;
my $plOutput;
my $minMAF;
my $maxMAF;
my $kCT;  #keepComplementaryTransversions;
my $license_opt;

$| = 1; # Reduce STDOUT buffer size

# Constants for annotation hash, recommended not to change 0-7 as these are hard coded in several places
use constant ANNOT_SNP    => 0; # SNP name
use constant ANNOT_CHR    => 1; # Chromosome
use constant ANNOT_ALA    => 2; # Allele A	
use constant ANNOT_GPO    => 3; # Sex averaged genetic position
use constant ANNOT_FRE    => 4; # Allele frequency of A
use constant ANNOT_HET    => 5; # Heterozygosity
use constant ANNOT_RSN    => 6; # rs name
use constant ANNOT_POS    => 7; # Physical position
use constant ANNOT_ALB    => 8; # Allele B	
use constant ANNOT_TBO    => 9; # Reference TOP/BOT strand
use constant ANNOT_TOTAL_COLS => 10 ; # Total columns of annotation hash. Should be +1 of above

use constant UNIVERSAL_POPCOL => 8; # The first column of the frequencies in the annotation file

GetOptions(
    "data=s"=>\$data,
    "annotDir=s"=>\$annotDir,
    "annotFile=s"=>\$annotFile,
    "crlmm"=>\$crlmm,
    "noX"=>\$noX,
    "randomSNP"=>\$randomSNP,
    "seed=s"=>\$seed,
    "pedfile=s"=>\$pedfile,
    "prog=s"=>\@prog,
    "whichSamplesFile=s"=>\$whichSamplesFile,
    "whichSamplesList=s"=>\$whichSamplesList,
    "callDir=s"=>\$callDir,
    "callFile=s"=>\$callFile,
    "chip=i"=>\$chip,
    "binsize:s"=>\$binsize,
    "outputDir=s"=>\$outputDir,
#    "keepME"=>\$keepME,
	"MendelErrors=s"=>\$MEO,
    "removeWFHBS:s"=>\$removeWFHBS,
    "removeAIS"=>\$removeAIS,
    "pop=s"=>\$pop,
    "popCol=i"=>\$popCol,
    "minDist:s"=>\$minDist,
    "freq"=>\$freq,
    "regions=s"=>\$regions,
    "regionsFile=s"=>\$regionsFile,
    "popHetTest:s"=>\$popHetTest,
    "bestPopTest"=>\$bestPopTest,
    "actg"=>\$actg,
    "fileRemoveSNPs=s"=>\$fileRemoveSNPs,
    "fileKeepSNPs=s"=>\$fileKeepSNPs,
    "help"=>\$help,
    "plSim=i"=>\$plSim,
    "plOutput:s"=>\$plOutput,
    "minMAF=f"=>\$minMAF,
    "maxMAF=f"=>\$maxMAF,
    "kCT"=>\$kCT,
    "keepComplementaryTransversions"=>\$kCT,
    "license"=>\$license_opt, 
) or print_usage("Error in parsing options: unknown option/s given (see top of this message for more details).");
if(@ARGV > 0) { print_usage("Unused parameters in command line: " . join("\t", @ARGV)) };

if(defined($help)){
    pod2usage(
        -verbose => 2, 
        -output => \*STDOUT, 
        -width => 50,
        -noperldoc => ! -t STDOUT,
    );
}


# Although there are ways to print out the perldoc directly, this requires higher versions of Perl to 
# operate properly, so just copy-paste the perldoc into a Here-Doc
print <<BANNER;

                    ------------------------------- LINKDATAGEN -------------------------------
                    now incorporating:
                    LINKDATAGEN_AFFY
                    LINKDATAGEN_ILLUMINA
                    LINKDATAGEN_MPS
                    M Bahlo, CJ Bromhead
        
                    with help from Tom Scerri, Katherine Smith, Rick Tankard and Luke Gandolfo.
                    ---------------------------------------------------------------------------

If you use LINKDATAGEN please acknowledge by citing:
Bahlo M, Bromhead CJ. Generating linkage mapping files from Affymetrix SNP chip data. Bioinformatics 2009;25(15):1961-2.

If you use LINKDATAGEN -data m (for MPS data) please acknowledge by citing:
Smith KR, Bromhead CJ, Hildebrand MS, Shearer AE, Lockhart PJ, Najmabadi H, Leventer RJ, McGillivray G, Amor DJ, Smith RJ, Bahlo M. 
Reducing the exome search space for Mendelian diseases using genetic linkage analysis of exome genotypes. Genome Biology 2011;12:R85.

This script is free software that is licensed under GPL v2 and comes with no warranty. For details, run the -license option. 

This script updated on $last_changed_date (v$script_revision).

BANNER

if(defined($license_opt)) {
    print $license;
    exit(1);
}

print "Started at ";
&print_time();
print "\n";


if(!defined($data) || $data !~ /^[aim]$/){
	print_usage("Must define a data type (a, i or m).");
}

if((defined($callDir) && defined($callFile))|| (!defined($callDir) && !defined($callFile))){
	print_usage("Must specify either a callDir (with -data i) or a callFile (with -data a or -data m).");
}

if($data eq "m" || $data eq "a"){
	if(!defined($callFile)){
		print_usage("This data type (-data m or -data a) requires the -callFile option (a single file containing all the genotype data).");
	}
	
	if(defined($whichSamplesList)){
		print_usage("This data type (-data m or -data a) requires the -whichSamplesFile option, not the whichSamplesList option.");
	}
}

#if($data eq "i"){
#	if(!defined($callDir)){
#		print_usage("This datatype (-data i) should use the -callDir option (directory specification for multiple files containing the genotype data).");
#	}
#}

if(!defined($MEO)){

	print_usage ("Must define action for Mendel errors:  -MendelErrors [ keep | removeSNPs | removeGPF ]");

} elsif (!($MEO eq "keep" || $MEO eq "removeSNPs" || $MEO eq "removeGPF")) {

	print_usage ("Aborted due to wrong value for mendelErrors. Only keep, removeSNPs, or removeGPF allowed");

}

if((!defined($whichSamplesFile) && !defined($whichSamplesList)) || (defined($whichSamplesFile) && defined($whichSamplesList))){
	print_usage("Must specify one, and only one, of either option -whichSamplesFile or -whichSamplesList file.");
}

if(!defined($pedfile)){
		print_usage("Must define pedfile.");
}

if( !(scalar(@prog) > 0  || defined($popHetTest) || defined($freq)) xor defined($bestPopTest) ) {
		print_usage("Must define EITHER \[ -prog and/or -popHetTest and/or -freq \] OR  \[ -bestPopTest \] alone.");
}

if (scalar(@prog) > 0) {

	@prog = split(/,/,join(',',@prog));

	foreach my $prog (@prog) {

		$progs{$prog} = 1;
	}
}

if (defined($actg)) {

	print_usage("Option -actg is not yet working.  De-select this option and try again.");

	if($data eq "a") {
		print_usage("You have specified both \"-data a\" and \"-actg\", however these options are currently mutually exclusive. Option \"-actg\" only works with \"-data i\" or \"-data m\".");
	}

	if($progs{"me"} == 0) {

#		print_usage("Option -actg has been specified, but this only works for merlin output which has not been selected. Can not proceed as I don't know what you want me to do.\n");
	}

#	print "Output of genotypes will be in acgt format\n";
} 

if(defined($chip) && defined($annotFile)) {
	print_usage("Must only define -chip or -annotFile, not both.");
}

if(defined($chip)){
	if(!defined($annotDir)){
		print_usage("If -chip is defined, then -annotDir must also be defined to specify the location of the annotation file; otherwise use -annotFile instead of -chip to explicitly name the annotation file.");
	}
}
elsif (!defined($annotFile)) {
	print_usage("Must define -chip or -annotFile.");
}

if (defined($binsize)) {
	if ($binsize !~ /^\d*\.?\d+$/ || $binsize eq "") {

		print_usage("Value for -binsize must be a number >= 0.0 whereas your value for -binsize was: $binsize");
	}
}

else {

	$binsize=0.3;  # size of bin from which to extract ONE marker, to be made into an option. used to reduce size of files.
}

if (defined($minDist)) {
	if ($minDist !~ /^\d*\.?\d+$/ || $minDist eq "") {

		print_usage("Value for -minDist must be a number >= 0.0 whereas your value for -minDist was: $minDist");
	}
}

else {

	if ($binsize * 0.5 > 0.2) {

		$minDist = 0.2;
	}

	else{

		$minDist = $binsize * 0.5;
	}
}

if (defined($minMAF) && ($minMAF < 0 || $minMAF > 0.5)) {

	print_usage("minMAF must be a number >= 0.0 and <= 0.5; your given value for minMAF was: $minMAF");
}

if (defined($maxMAF) && ($maxMAF < 0 || $maxMAF > 0.5)) {

	print_usage("maxMAF must be a number >= 0.0 and <= 0.5; your given value for maxMAF was: $maxMAF");
}

if ((defined($minMAF) && defined($maxMAF)) && ($minMAF > $maxMAF)) {

	print_usage("minMAF must be less than or equal to maxMAF; your given value for minMAF was $minMAF and for maxMAF was $maxMAF");
}
#printf ("Markers will be forced to be at least %1.2f apart.\n",  $minDist);

if (defined($plOutput) && !($plOutput eq "allSamples" || $plOutput eq "makeAllFounders")) {

	print_usage("plOutput can only take values \"allSamples\" or \"makeAllFounders\"; your value for plOutput was: $plOutput");
}

my %regions;

if(defined($regions || $regionsFile)){

	my $tempRegion;
	my @tempRegions;
	my @tempRegion;

	if(defined($regions)) {

		#note - command-line regions should be comma-separated and in the format "chr#:#####-#####", where:
		#the first # must be a number 1 to 22, or else the letter X, Y or M.
		#the remaining #### are numbers representing start and stop positions (bp).

		@tempRegions = split(/,/,$regions); #split the comma-separated command line regions.

		foreach $tempRegion (@tempRegions) {

			chomp($tempRegion);

			if ($tempRegion =~ /^(?:chr)?([1-9]|1[0-9]|2[0-2]|[XYM]):(\d+)-(\d+)$/ ) {

				push (@{$regions{$1}}, [$2, $3] );		# %regions{chr} references an array "of arrays of start and stop positions".
			}

			elsif ($tempRegion =~ /^(?:chr)?([1-9]|1[0-9]|2[0-2]|[XYM])$/) {

				print "\nWhole of chromosome $1 selected.";	#whole chromosome selected
				push (@{$regions{$1}}, [0, 999999999] );	# %regions{chr} references an array "of arrays of start and stop positions".
			}

			else {
				print_usage("Unrecognised format of region \"$tempRegion\" from command-line region(s): $regions");
			}
		}
	}

	if (defined ($regionsFile))	{

		copen(*REG,"<$regionsFile") or print_usage "Can't open region file: $regionsFile. $!"; # region file

		foreach $tempRegion (<REG>) {

			chomp($tempRegion);

			if ($tempRegion =~ /^(?:chr)?([1-9]|1[0-9]|2[0-2]|[XYM])(?:\s+|:)(\d+)(?:\s+|-)(\d+)(?:\s*)$/ ) {

                                push (@{$regions{$1}}, [$2, $3] );              # %regions{chr} references an array "of arrays of start and stop positions".
                        }

			elsif ($tempRegion =~ /^(?:chr)?([1-9]|1[0-9]|2[0-2]|[XYM])(?:\s*)$/) {

				print "\nWhole of chromosome $1 selected.";     #whole chromosome selected
				push (@{$regions{$1}}, [0, 999999999] );        # %regions{chr} references an array "of arrays of start and stop positions".
			}

			else {
				print_usage("Unrecognised format of region \"$tempRegion\" from regionsFile: $regionsFile");
			}
		}
	}

	print "\n\nThe following regions to output have been selected from command line and/or file:\n\n";

	foreach my $tempChr (%regions) {

		foreach $tempRegion (@{$regions{$tempChr}}) {

			print "Chr $tempChr Start $tempRegion->[0] End $tempRegion->[1]\n";
		}
	}
}

if(defined($seed)) {

	if ($seed !~ /^-?\d+$/) {

		print_usage("Value for -seed must be an integer.");
	}

	elsif (!defined($randomSNP)) {

		print_usage("If using -seed then you must use \"-randomSNP\" as well, otherwise remove -seed from your command line.");
	}
}

elsif (defined($randomSNP)) {

	$seed = 12345;
}

if(!defined($randomSNP)){

#	print "\nThe SNP with greatest heterozygosity will be chosen from each interval.\n";
}

else {

	srand($seed);
#	print "\nA random SNP will be chosen in each interval with seed \"$seed\".  This seed can be changed with the -seed { integer } option.\n";
}

if(defined($removeWFHBS)){

	if ($removeWFHBS eq "u" || $removeWFHBS eq "") {

		$removeWFHBS = "u";
#		print "\nThe union of \"within-family\ homozygosity-by-state\" SNP markers across families will be removed.\n";
	}

	elsif ($removeWFHBS eq "i") {

#		print "\nThe intersection of \"within-family\ homozygosity-by-state\" SNP markers across families will be removed.\n";
	}

	else {
		print_usage("Valid options for -removeWFHBS are only \"u\" (union; the default) or \"i\" (intersection) for removal of \"within-family homozygosity-by-state\" SNP markers across families. Your option was: $removeWFHBS");
	}
}
else {

#	print "\n\"Within-family homozygosity-by-state\" SNP markers will NOT be removed.\n";
}

if (defined($popHetTest)) {

	if ( ! ($popHetTest eq "" || $popHetTest eq "summary" || $popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") ) {

		print_usage("Valid options for -popHetTest are only \"summary\" (default), \"verbose\", \"perChr\" or \"perChrVerbose\". Your option was: $popHetTest");
	}
	
	if ($popHetTest eq "") {

		$popHetTest = "summary";
	}
}

#  -----------Annotation files -------------------
my $chipfile;
my $chipfreqs;
if($data eq "a"){
    my @chips=("50KXba","50KHind","250KSty","250KNsp","5","6","500KNspSty");
	if(defined($annotFile)) {
		if(defined($annotDir)){
			$chipfile = $annotDir."/".$annotFile;
		}
		else{
			$chipfile = $annotFile;
		}
	}
	else{
		if($chip >= 1 && $chip <= 7){
			$chipfile=$annotDir."/annot".$chips[$chip-1].".txt";
		}
		else{
			print_usage("If -data a is specified, then -chip must be 1, 2, 3, 4, 5, 6 or 7.");
		}
	}
}
if($data eq "i"){
    my @chips=("370","610","660","Cyto12","OmniEx","1m","HapMap2U","HapMap3U");
	if(defined($annotFile)) {
		if(defined($annotDir)){
			$chipfile = $annotDir."/".$annotFile;
		}
		else{
			$chipfile = $annotFile;
		}
        # determine annotation type of Illumina data
	}
	else{
		if($chip >= 1 && $chip <= 8){
			$chipfile=$annotDir."/annot".$chips[$chip-1].".txt";
            $annotFileIlluminaType = 'old_illumina';
		}
		else{
			print_usage("If -data i is specified, then -chip must be 1, 2, 3, 4, 5, 6, 7 or 8.\n");
		}
	}
    # Check type of Illumina annotation file
    copen(*IN2,"<$chipfile") or print_usage "Can't open annotation file $chipfile. $!";
    while(<IN2>) {
        if(/^Probe_set_ID/) {
            # This is the old annotation type
            $annotFileIlluminaType = 'old_illumina';
            last;
        }
        if(/Minimum LINKDATAGEN revision/) {
            # This is the new annotation type
            last;
        }
        unless(/^#/) {
            print_usage("Input annotation file not recognised as valid, being neither the June 2011 or newer formats.");
        }
    }
    close(IN2);
}
if($data eq "m"){
    my @chips=("HapMap2U","HapMap3U");
	if(defined($annotFile)) {
		if(defined($annotDir)){
			$chipfile = $annotDir."/".$annotFile;
		}
		else{
			$chipfile = $annotFile;
		}
	}
	else{
		if($chip >= 1 && $chip <= 2){
			$tmp=$chips[$chip-1];
			$chipfile=$annotDir."/annot".$tmp.".txt";
		}
		else{
			print_usage("If -data m is specified, then -chip must be 1 or 2.\n");
		}
	}
}



# -------------- populations -------------------
my %pops=(
     CEU => 1,
     ASW => 2,
     CHB => 3,
     CHD => 4,
     GIH => 5,
     JPT => 6,
     LWK => 7,
     MEX => 8,
     MKK => 9,
     TSI => 10,
     YRI => 11,
);
my @hapmap2 = qw (CEU CHB JPT YRI); 

my $key;

print "\n";

#foreach $key (keys %pops){
#	print "(",$key,",",$pops{$key},")\n";
#}

if((defined($popCol) + defined($pop) + defined($bestPopTest)) > 1){

		print_usage("Options -popCol, -pop and -bestPopTest are mutually exclusive, they can not be specified together.  Please specify only one.");
}

if(defined($popCol)){
	
	--$popCol;
	print "User (USR) defined population allele frequencies have been selected from column $popCol (zero-based; i.e. 1 less than your given value)\n";
	$pop = "USR";
}
elsif(!defined($pop)){
	#set pop to the default of CEU (important note - this is the first population in the annotation file which matters due to the hard-wired nature of LINKDATAGEN)
    if($annotFileIlluminaType eq 'old_illumina') {
        $popCol = 7; 
	} elsif($data eq "m" || $data eq "i"){
	    $pop = "CEU";
        $popCol = UNIVERSAL_POPCOL;
	} elsif($data eq "a"){
		$popCol=6;
	}

	if (!defined($bestPopTest)) {
	
		print "No population has been specified for the allele frequencies: Caucasians (CEU) chosen by default\n";
	}

#	$pop = "CEU";
}
else{
	print "Population $pop has been specified for the allele frequencies.\n";
    unless(defined($pops{$pop})) {
        print_usage "The population $pop is not known by LINKDATAGEN";
    }
	if($data eq "m"){
		if((defined($chip) && $chip==1)|| (defined($annotFile) && $annotFile=~/annotHapMap2/)){
			foreach my $hapm (@hapmap2){
				if($pop eq $hapm){ # is one of the hapmap2 pops and ok
					goto POPCOL;
				}
			}
			print_usage("For MPS data, currently the only choice is the CEU, CHB, JPT or YRI for HapMap2 data.");
		}
		if((defined($chip) && $chip==2)||(defined($annotFile) && $annotFile=~/annotHapMap3/)){
			foreach $key (@hapmap2){
				if($pop eq $key){ # is one of the hapmap2 pops
					print "WARNING: Use annotHapmap2.txt file instead, for a much greater selection of SNPs.\n";
				}
			}
		}
	}
    POPCOL: 
    if($annotFileIlluminaType eq 'old_illumina') {
        $popCol= 6 + $pops{$pop};
    } elsif($data eq "i" || $data eq "m"){
		$popCol= UNIVERSAL_POPCOL - 1 + $pops{$pop};
	} elsif($data eq "a"){
		$popCol= 5 + $pops{$pop};
	}
}
# chromosome numbers
$numchrom = 23;        #CATHEDIT: Fix for if you want to exclude X chromosome
if(defined($noX)) {
	$numchrom = 22; 
}

# ---------OUTPUTDIRS -----------------
my $merlinDir;
my $allegroDir;
my $prestDir;
my $morganDir;
my $plinkDir;
my $cpDir;
my $beagleDir;
my $festimDir;
my $relateDir;
my $tablesDir;
my $fastPHASEDir;

if($progs{"all"} == 1){
	%progs = ("al"=>1,"fp"=>1,"me"=>1,"pr"=>1,"cp"=>1,"pl"=>1,"mo"=>1,"be"=>1,"fe"=>1,"re"=>1,"all"=>1);
}

if(defined($outputDir)){
	$outputDir = $outputDir . "_";
}
else {
	$outputDir = "";
}

$tablesDir = $outputDir . "tables";	#this will be the directory to store all the tables, e.g. QC, Mendel errors etc...
mkdirnice($tablesDir);

if($progs{"al"} == 1){
	$allegroDir=$outputDir . "allegro";
	print "The output dir for the ALLEGRO files will be $allegroDir/\n";
	mkdirnice($allegroDir);
}
if($progs{"me"} == 1){
	$merlinDir=$outputDir . "merlin";
	print "The output dir for the MERLIN files will be $merlinDir/\n";
	mkdirnice($merlinDir);
	mkdirnice($merlinDir. "/genome");
}
if($progs{"pr"} == 1){
	$prestDir=$outputDir . "prest";
	print "The output dir for the PREST files will be $prestDir/\n";
	mkdirnice($prestDir);
}
if($progs{"mo"} == 1){
	$morganDir=$outputDir . "morgan";
	print "The output dir for the MORGAN files will be $morganDir/\n";
	mkdirnice($morganDir);
}
if($progs{"pl"} == 1){
	$plinkDir=$outputDir . "plink";
	print "The output dir for the PLINK files will be $plinkDir/\n";
	mkdirnice($plinkDir);
}
if($progs{"cp"} == 1){
	$cpDir=$outputDir . "complete";
	print "The output dir for the \"complete\" format files will be $cpDir/\n";
	mkdirnice($cpDir);
	}
if($progs{"be"} == 1){
	$beagleDir=$outputDir . "beagle";
	print "The output dir for the BEAGLE files will be $beagleDir/\n";
	mkdirnice($beagleDir);
}
if($progs{"fe"} == 1){
	$festimDir=$outputDir . "festim";
	print "The output dir for the FESTIM files will be $festimDir/\n";
	mkdirnice($festimDir);
}
if($progs{"re"} == 1){
	$relateDir=$outputDir . "relate";
	print "The output dir for the RELATE files will be $relateDir/\n";
	mkdirnice($relateDir);
}
if($progs{"fp"} == 1){
	$fastPHASEDir=$outputDir . "fastPHASE";
	print "The output dir for the fastPHASE files will be $fastPHASEDir/\n";
	mkdirnice($fastPHASEDir);
}

# ------------BRLMM/CRLMM genotype calls -----------
if($data eq "a"){
	if(!defined($crlmm)){
		print "BRLMM genotype calls. Genotypes are one of {-1,0,1,2}.\n\n";
		$brlmm=0;
	}
	else{
		print "CRLMM genotype calls. Genotypes are one of {0,1,2,3}.\n\n";
		undef($brlmm);
	}
}
if($data eq "i" || $data eq "m"){
	$brlmm=0;
}


# ----------frequencies --------------
if(defined($freq)){
	print "Allele frequencies estimated from the founders in the pedigree will be written to alleleFreqs.txt\n";
}

# -----global vars------
my %genos_orig=(); # hash of arrays with hash containing ALL SNPs and array containing genotypes of people
my %genosOrigArchive=();
my %annot_orig=(); # hash of arrays with hash containting ALL SNP ids and array  containing annotation data
my %genos=(); # hash of arrays with hash containing chosen subset of SNPs and array containing genotypes of people
my %hash_header_index=();  #index of header with which to recognise column names
#my %annot=(); # hash of arrays with hash containting chosen subset of SNP ids and array  containing annotation data
my @ped; # pedigree array
#my @whichsamples=(0,0,0,0,0,0,0,0,14,12,11,13); # these are the samples that need to be put in the .pre file, others are going to be missing
my @whichsamples=();
my @whichsamples_list=(); #CATHEDIT
my $pedno=0; # no of individuals in the pedigree
my @chr_snp_cnt_orig=(); # numbers of Snps on each chromosome for ALL Snps
my @chr_order_snp_orig=(); # array (chr) of array (genetic map) of array (all annot data) of snps ordered by 1) chromosome and 2) genetic map pos
my @chr_snp_cnt=(); # numbers of Snps on each chromosome for subset of snps
my @chr_order_snp=(); # array (chr) of array (genetic map) of array (all annot data) of snps ordered by 1) chromosome and 2) genetic map pos
my @tempMB;
my $pfile;
my $command;
my %bin_het_hash=();
my @which=();
my @whichlist=(); #CATHEDIT
my @whichnames=();


@tempMB=split(/\./,$pedfile);
$pfile=$tempMB[0]; # pedigree name

my $line;
my @temp;
my %temp_chr_hash=();

# order of the opening of these files is IMPORTANT!!!!!
#read in the whichsamples file
my @arr4; #for whichSamplesFile
my @arr5;	#for whichSamplesList
if(defined($whichSamplesFile)) {
	copen(*IN4,"<$whichSamplesFile") or print_usage "Can't open whichSamplesFile $whichSamplesFile file. $!";
	@arr4=<IN4>;
	print "\nReading in the whichsamples file...\n\n";
	read_in_whichsamples();
	undef(@arr4);
	close(IN4);
}
else{ #CATHEDIT: Get @whichsamples() and @whichsamples_list() instead
	copen(*IN5, "<$whichSamplesList") or print_usage "Can't open whichsamples list $whichSamplesList file. $!";
	@arr5=<IN5>;
	print "\nReading in the whichSamplesList file...\n\n";
	read_in_whichSamplesList();
	undef(@arr5);
	close(IN5);
}

print "\nReading in the pedigree data ...\n";

copen(*IN3,"<$pedfile") or print_usage "Can't open pedigree file $pedfile. $!"; # pedigree file
my @arr3;
@arr3=<IN3>;
read_in_ped();
undef(@arr3);
close(IN3);


print "\nMatching pedigree data to whichsample data ...\n";

for(my $i=0;$i<=$#whichsamples;$i++){
	if($whichsamples[$i]!=0){
		push(@whichnames, $ped[$i][0]);
		push(@whichnames,$ped[$i][1]); 
		push(@which,$whichsamples[$i]); # which gives the nth sample which 				
		if(defined($whichSamplesList)){
			push(@whichlist, $whichsamples_list[$i]);
		}
		## +1 because SNP is also in this line so counting from 1 to no_samples in whichsamples is fine and does not need to be adjusted
	} #  which is of size <= whichsampels as all missing inds are not included
}

print "\nSorting pedigree data into families...\n";

my %fam=();		#hash of all unique families, each pointing to a hash of all samples within each family defined by the pedigree (regardless of genotyping)
my %fam_i=();	#hash of all unique families, each pointing to an array of only genotyped samples with their unique position in the overall genotype matrix
my $familiesTotal;
my %famgenocnt;	#number of individuals genotyped in this family
my $family;
my $individual;
my $i;

for($i=0;$i<=$#whichsamples;$i++){
		$fam{$ped[$i][0]}{$ped[$i][1]}=$i;
}

$familiesTotal = scalar keys (%fam);

print "\nTotal number of families = $familiesTotal\n";

foreach $family (keys %fam){

	foreach $individual (keys %{$fam{$family}}){

		if($whichsamples[$fam{$family}{$individual}]!=0){
		
			$famgenocnt{$family}+=1;	#add 1 to the number of genotyped individuals in this family
			
			for($i=0;$i<=$#which;$i++){
			
				if($whichsamples[$fam{$family}{$individual}]==$which[$i]){

					push(@{$fam_i{$family}},$i);
					last;
				}
			}
			
#			print keys @{$fam_i{$family}};
		}
	}
}


############

#New code to loop over the whole thing for the popTest.
#

my %bestPopTest;
my %autosomeCountSNPs;
my %popAutosomeMeanHet;
my $bestPop;
my $popColOriginal = $popCol;
my @pops = ("CEU", "ASW", "CHB", "CHD", "GIH", "JPT", "LWK", "MEX", "MKK", "TSI", "YRI");
my $lastPop = $popCol + $#pops;

if (defined($bestPopTest)) {

	$popHetTest = "summary";

    # for $pop (@pops) {
	for (; $popCol <= $lastPop; $popCol++) {
		works();
	}

	open(BPT, ">$tablesDir" . "/bestPopTest.txt") || die print "Can't open bestPopTest.txt\n";

	print "\nDisplaying the results of -bestPopTest in the table below.  Shown are the chi-sq values derived from the autosome (not genome) only:\n\n";
	print "\nFamily\tMember\tCEU\tASW\tCHB\tCHD\tGIH\tJPT\tLWK\tMEX\tMKK\tTSI\tYRI\t\tBEST";
	print "\n-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t\t-------";

	print BPT "Displaying the results of -bestPopTest in the table below.  Shown are the chi-sq values derived from the autosome (not genome) only:\n\n";
	print BPT "\nFamily\tMember\tCEU\tASW\tCHB\tCHD\tGIH\tJPT\tLWK\tMEX\tMKK\tTSI\tYRI\t\tBEST";
	print BPT "\n-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t\t-------";

	for (my $i=0; $i<$pedno; $i++){

	my $bestPop = $popColOriginal;

		for(my $j=0;$j<=$#which;$j++){ # genotyped individual

			if($which[$j]==$whichsamples[$i]){

				printf ("\n%s\t%s", substr($ped[$i][0], 0, 7), substr($ped[$i][1], 0, 7));
				printf BPT ("\n%s\t%s", $ped[$i][0], $ped[$i][1]);

				for ($popCol = $popColOriginal; $popCol <= $lastPop; $popCol++) {

					if (defined(@{$bestPopTest{$popCol}}[$i])) {

						printf ("\t%1.3f", @{$bestPopTest{$popCol}}[$i]);
						printf BPT ("\t%1.3f", @{$bestPopTest{$popCol}}[$i]);
						
						if (@{$bestPopTest{$popCol}}[$i] < @{$bestPopTest{$bestPop}}[$i]) {

							$bestPop = $popCol;
						}
					}
					else {
						print "\tN/A";
						print BPT "\tN/A";
					}
				}

				print "\t\t$pops[$bestPop - $popColOriginal]";
				print BPT "\t\t$pops[$bestPop - $popColOriginal]";

				if ( length($ped[$i][0]) > 7) {

					print "\tFamily name truncated.";
				}
				
				if ( length($ped[$i][0]) > 7) {
				
					print "\tMember's name truncated.";
				}
				
				if ( length($ped[$i][0]) > 7 || length($ped[$i][1]) > 7) {

					print "\tSee the tab-delimited popHetTest.txt table for full family/member names.";
				}
			}
		}
	}

	print "\n-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t\t-------";
	print "\nNumber of SNPs";
	print BPT "\n-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t\t-------";
	print BPT "\nNumber of SNPs\t";


	for ($popCol = $popColOriginal; $popCol <= $lastPop; $popCol++) {

		if (defined($autosomeCountSNPs{$popCol})) {

			print "\t$autosomeCountSNPs{$popCol}";
			print BPT "\t$autosomeCountSNPs{$popCol}";
		}
		else {

			print "\t0";
			print BPT "\t0";
		}
	}

	print "\nMean heterozyg";
	print BPT "\nMean heterozyg\t";

	for ($popCol = $popColOriginal; $popCol <= $lastPop; $popCol++) {

		if (defined($popAutosomeMeanHet{$popCol})) {

			printf ("\t%1.3f", $popAutosomeMeanHet{$popCol});
			printf BPT ("\t%1.3f", $popAutosomeMeanHet{$popCol});
		}
		else {

			print "\tN/A";
			print BPT "\tN/A";
		}
	}

	print "\n-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t";
	print BPT "\n-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t-------\t";

	if (defined($removeWFHBS)) {
	
		print "\nWARNING: \"Within-family homozygosity-by-state\" SNP markers have been removed.  This will make the available SNPs more heterozygous than by chance alone, thus adversely effecting this analysis.  It is recommended NOT to use the \"-removeWFHBS\" option when running this analysis.\n";
		print BPT "\nWARNING: \"Within-family homozygosity-by-state\" SNP markers have been removed.  This will make the available SNPs more heterozygous than by chance alone, thus adversely effecting this analysis.  It is recommended NOT to use the \"-removeWFHBS\" option when running this analysis.\n";
	}

	close (BPT);

	print "\n\n\tASW : African ancestry in Southwest USA\n";
	print "\tCEU : Utah residents with Northern and Western European ancestry from the CEPH collection\n";
	print "\tCHB : Han Chinese in Beijing, China\n";
	print "\tCHD : Chinese in Metropolitan Denver, Colorado\n";
	print "\tGIH : Gujarati Indians in Houston, Texas\n";
	print "\tJPT : Japanese in Tokyo, Japan\n";
	print "\tLWK : Luhya in Webuye, Kenya\n";
	print "\tMEX : Mexican ancestry in Los Angeles, California\n";
	print "\tMKK : Maasai in Kinyawa, Kenya\n";
	print "\tTSI : Toscans in Italy\n";
	print "\tYRI : Yoruba in Ibadan, Nigeria (West Africa)\n";

	print "\nPlease consult the above table and decide upon the best population to use for your sample, and then re-run linkdatagen using the -pop option.\n\n";

	config();

	print "\nFinished at ";
	&print_time();

	exit();
}
else {

	works();
}

############


sub works {

my $line_cnt;

undef(%annot_orig);
undef(@chr_snp_cnt_orig);

print "\nOpening annotation file $chipfile ...\n";

copen(*IN2,"<$chipfile") or print_usage "Can't open annotation file $chipfile. $!";# annotation file with cM, all freqs, hets etc

#if ($pop ne "USR") {
#	$pop = $pops[$popCol - $popColOriginal];
#}

#print "Reading in the genome annotation data using population allele frequencies from column $pop ...\n";

if (defined($bestPopTest)) {
	print "Reading in the genome annotation data using population allele frequencies from column $pops[$popCol - $popColOriginal] ...\n";
} elsif (defined($pop)) {
	print "Reading in the genome annotation data using population allele frequencies from column $pop ...\n";
#} elsif (defined($pop)) {
#	print "Reading in the genome annotation data using population allele frequencies from column " . ($pop eq "USR" ? $popCol + 1 : $pop ) . " ...\n";
} else {
	print "Reading in the genome annotation data using CEU population allele frequencies (default selection)...\n";
}



if (defined($popHetTest)) {

	print "Skipping SNPs with a minor allele frequency (MAF) < 0.4 as -bestPopTest (or -popHetTest) has been selected.  This is important for the subsequent chi-square test(s) ...\n";
}

if($data eq "a"){
	$line_cnt = read_in_annot_affy();
}
if($data eq "i"){
    if($annotFileIlluminaType eq 'old_illumina') {
        $line_cnt = read_in_annot_illumina_legacy();
    } else { 
        $line_cnt = read_in_annot_mps_illuminaSNP();
    } 
}
if($data eq "m"){
	$line_cnt = read_in_annot_mps_illuminaSNP();
}


if (defined($fileRemoveSNPs)) {

	print "Reading fileRemoveSNPs $fileRemoveSNPs...\n"; 

	copen(*FRS,"<$fileRemoveSNPs") or print_usage "Can't open $fileRemoveSNPs. $!";
	
	removeSNPs();

	close(FRS);
}

if (defined($fileKeepSNPs)) {

	print "Reading fileKeepSNPs $fileKeepSNPs...\n"; 

	copen(*FKS,"<$fileKeepSNPs") or print_usage "Can't open $fileKeepSNPs. $!";
	
	keepSNPs();

	close(FKS);
}


if ($line_cnt == 0) {

	if (!defined($bestPopTest)) {

		print_usage("Can not continue as the annotation file contains no usable SNPs for the selected population.");
	}
	else {
		return;
	}
}

close(IN2);

if (! (%genosOrigArchive) ) {

	# input data manipulation
	print "\nReading in the genotyping data ...\n\n";
	read_in_brlmm();

	print "There are " . scalar (keys %genos_orig) . " keys.\n";
}

else {

	undef(%genos_orig);

	print "Retrieving genotype data from memory.\n";
	%genos_orig = %genosOrigArchive;
	print "There are now " . scalar (keys %genos_orig) . " keys.\n";
}


# verify and fix data so that all SNPs have annotation data and genotype data
check_geno_markers();

#print_data_test();

# need to recode all the genotype data .. into brlmm data
if($data eq "i"){
	recode_illumina_data();
}

# test print data function
#print_data_test();

# error checks - mainly for Allegro & MORGAN, merlin knocks them out itself

print "\nChecking genotyping data for Mendelian errors ...\n";

mendelian_errors();

# remove all "within-family homozygosity-by-state" SNP marker.  Maybe intersection or union across families.
if(defined($removeWFHBS)){

	print "\nChecking genotype data for \"within-family homozygosity-by-state\" SNP markers ...\n";

	removeWFHBS();
}

print "\nSorting SNPs by chromosome ...\n\n";

split_snps_to_chromos();
print_ordered_snps();

# printing out allele frequencies for all 
if(defined($freq)){
	allele_freqs();
}

# functions to reduce datafiles to fewer markers

if($binsize==0){

	print "\nBinsize=0, so all SNPs have been selected...\n";
}

else {

	if (defined ($randomSNP)) {

		print "\nReducing dataset by selecting a random SNP from each bin ...\n";
	}

	else {

		print "\nReducing dataset by selecting the SNP with highest heterozygosity from each bin ...\n";
	}
}

open(OUT1, ">$tablesDir" . "/missingMarkers.txt") || die print "Can't open logfile with missingMarkers.txt\n";

select_per_bin();

close(OUT1);

print "\nPrinting selected SNPs to selectedSNPs.txt ...\n";

open(SELECTEDSNPS, ">$tablesDir/selectedSNPs.txt") or print_usage("Cannot open file $tablesDir/selectedSNPs.txt, $!");

foreach $key (keys %genos_orig) {

	print SELECTEDSNPS $key, "\n";
}
	
close(SELECTEDSNPS);

if (defined($popHetTest)) {

	$bestPopTest{$popCol} = popHetTest();
}

# test ---------
#print_genos_test();

undef(%genos_orig);
undef(%annot_orig);


###############
#

}

#end of new loop best pop test

#
###############


# output files
if($progs{"all"} == 1) {
	print "Selected ALL programs\n";
}

#merlin files
if($progs{"me"} == 1){
	print "Printing out merlin .pre files..\n";
	print_pre_merlin();
	print "Printing out merlin .map files...\n";
	print_map_merlin();
	print "Printing out merlin .freq files....\n";
	print_freq_merlin();
	print "Printing out merlin .dat files...\n";
	print_dat_merlin();
	print "Printing out merlin .in files with which to run each chromosome...\n";
	print_merlin_run_files();
}

# allegro files
if($progs{"al"} == 1){
	print "Printing out allegro .pre files..\n";
	print_pre_allegro();
	print "Printing out allegro .dat files..\n";
	print_dat_allegro();
	print "Printing out allegro .in files..\n";
	print_allegro_run_files();
}

# prest files
if($progs{"pr"} == 1){
	print "Printing .idx files for PREST\n";
	print_PREST_idx();
}

# Morgan files
if($progs{"mo"} == 1){
	print "Printing .par, .map and .ped files for MORGAN\n";
	print_MORGAN();
}

if($progs{"pl"} == 1){
	print "Printing out .map and .ped file for PLINK\n";
	print_PLINK();
}

if($progs{"cp"} == 1){
	print "Printing out the genotypes in complete style\n";
	print_CP();
}

if($progs{"be"} == 1){
	print "Making Beagle files\n";
	print_BEAGLE(); 
}
if($progs{"fe"} == 1){
	print "Printing data and map files for FEstim\n";
	print_FESTIM(); 
}

if($progs{"re"} == 1){
	print "Printing files for RELATE\n";
	print_RELATE();
}
if($progs{"fp"} == 1){
	print "Printing files for fastPHASE\n";
	print_fastPHASE();
}


print "\nAdditional output files are:\n";
print "missingMarkers.txt - contains bins with no SNP representation.\n";
print "chrX_SNPs.txt - contains SNPs on the X chromosome.\n";
print "chrY_SNPs.txt - contains SNPs on the Y chromosome.\n";
print "chrMT_SNPs.txt - contains SNPs on the mitochondrial genome.\n";
print "mendelErrors.txt - contains SNPs with Mendelian errors ordered by chromosome and map positions. Should highlight gross cytological abnormalities.\n";


config();

print "\nFinished at ";
&print_time();


# ---------sub routines--------------
sub print_usage{
	die("\nLINKDATAGEN has aborted.\n\n$_[0]\n\nUse -help for further instructions, or refer to the LINKDATAGEN manual.\n");
}

sub print_time() {
	my @time = localtime(time);
	if(length($time[1]) < 2) {
		$time[1] = "0".$time[1];
	}
	if(length($time[0]) < 2) {  #tim tells me the time!
		$time[0] = "0".$time[0];
	}
	my $tim = $time[2].":".$time[1].":".$time[0]."\n";
	print $tim;
}

sub copen {
    # Compressed open, supports opening of gzip and bzip2 compressed files
    # Supported usage, can be used like normal open (with * or scalar) in these forms only:
    #   copen *FILEHANDLE,EXPR
    #   copen *FILEHANDLE,MODE,EXPR
    #   copen $filehandle,MODE
    #   copen $filehandle,MODE,EXPR
    #my $fh = $_[0];
    my $mode;
    my $file;
    if(@_ == 2) {
        # copen FILEHANDLE,EXPR
        $_[1] =~ /^([<>])(.+)$/;
        $mode = $1;
        $file = $2;
    } elsif (@_ == 3) {
        # copen FILEHANDLE,MODE,EXPR
        $mode = $_[1];
        $file = $_[2];
    } else {
        die "copen function incorrectly used. (Program bug, please contact developer)"
    }
    # Check that we are not trying to open a .tar file
    if($file =~ /.\.tar(\.\w+)?$/ || $file =~ /\.t(b(z|z2|2)|gz)$/ ) {
        print_usage "Tar archive reading/writing is not supported. File: $file"; 
    }
    # change mode required for gzip and bzip2 files
    if($file =~ /\.gz$/) {
        # gzip file open
        unless($perliogzip) { print_usage "Perl module PerlIO::gzip not installed. Required for reading/writing gzip compressed files.\nOpen was attempted on: $file" }
        $mode = "$mode:gzip"; 
    } elsif ($file =~ /\.bz2$/) {
        unless($perlioviabzip2) { print_usage "Perl module PerlIO::via::Bzip2 not installed. Required for reading/writing bzip2 compressed files.\nOpen was attempted on: $file" }
        $mode = "$mode:via(Bzip2)";
    } 
    return( open($_[0], $mode, $file) );
}

sub mkdirnice {
    # Creates directories with easier to understand messages and better 
    # response to them not being created. 
    unless(mkdir $_[0]) {
        if($! eq 'File exists') {
            print "\nTried to create directory '$_[0]' but it already exists.\nOutput files from previous runs in this directory may be overwritten.\n\n";
        } else {
            print_usage "Directory '$_[0]' could not be created: $!.";
        }
    }
}

sub array_index{   #not sure about this sub routine, whether "return" is the right thing to say
	my $c;
	my %index_of_array1=();
	for($c = 0; $c < @_; $c++) {    #here @_ is the array we wish to index, and $_ is used to refer directly
		$index_of_array1{$_[$c]} = $c;			#to elements of the array
	}
	return %index_of_array1;
}

sub heterozygosity_from_af($) {
    # heterozygosity from allle frequency
    my $af = shift;
    if ($af < 0 || $af > 1) {
        print_usage "(Internal error) Allele frequency out of 0 to 1 range: $af. Possible error with annotation file.";
    }
    2 * $af * (1 - $af);
}

sub read_in_whichsamples{
	@whichsamples=split(/\s+/,$arr4[0]);
	print "whichsamples are: @whichsamples\n";
}

sub read_in_whichSamplesList{  #CATHEDIT: New subroutine to read in whichsamples list file.
	my $i=0;
	my @tmp;
	my $x = 1;
	for($i = 0; $i <=$#arr5; $i++) {
		@tmp = split(/\n/,$arr5[$i]);
		$whichsamples_list[$i] = $tmp[0];
		if($tmp[0] eq "0") {
			push(@whichsamples,0); 
		}
		else{
			push(@whichsamples,$x); 
			$x++;
		}
	}
	print "whichSamplesList\twhichsamples\n";
	for($i = 0; $i <=$#whichsamples_list; $i++) {
		print "$whichsamples_list[$i]\t$whichsamples[$i]\n";
	}
}


# There is a lot of header space for the illumina genotyping files.
# ditch the first 19 lines.
#Sample ID	SNP Name	SNP Index	Allele1 - Top	Allele2 - Top	GT Score	Theta	R	X	Y	X Raw	Y Raw	Log R Ratio	B Allele Freq
#
#1-EH	200003	1	A	A	0.8931	0.020	1.105	1.072	0.034	13377	1666	-0.0145	0.0000
#
sub read_in_brlmm{
#@which=();
my $i;
my $line_cnt=0;
my $k;
my @arr=();
my $t;
my $filename;
my @header=();
my @headarray=();
my $headerline;
my $count;
my $SNP_col;
my $alleleA_col;
my $alleleB_col;
my $sampleName_col;
my $sampleName;
my @num_rs_snps=();
my @hom_cnt=();
my @het_cnt=();
my $homrate;
my $hetrate;
my $missrate;
my @tmp=();
#my @whichnames=();
# affy vars 
my @nocall_cnt=();
my $str="SNP_A";
my @total_line_cnt=();
my $entry;
my $sum=0;
my $den;
#	for($i=0;$i<=$#whichsamples;$i++){
#		if($whichsamples[$i]!=0){
#			push(@whichnames, $ped[$i][0]);
#			push(@whichnames,$ped[$i][1]); 
#			push(@which,$whichsamples[$i]); # which gives the nth sample which 				
#			if(defined($whichSamplesList)){
#				push(@whichlist, $whichsamples_list[$i]);
#			}
#			## +1 because SNP is also in this line so counting from 1 to no_samples in whichsamples is fine and does not need to be adjusted
#		} #  which is of size <= whichsampels as all missing inds are not included
#	}
	for($k=0;$k<=$#which;$k++){
		$nocall_cnt[$k] = 0;
		$hom_cnt[$k] = 0;
		$het_cnt[$k] = 0;
	}

	if($data eq "i") {
	
		if (!defined($callFile)) {
		#my %trans=(AA,AC,AG,AT,CC,CG,CT,GG,GT,TT);
		# for illumina the whichsamples refer to the order of the files 
		for ($k=0;$k<=($#which);$k++) {
			$line_cnt=0;
			# need to find a prefix input setup for Illumina
			if(defined($whichSamplesFile)) {
				$filename=$callDir."_FinalReport".$which[$k].".txt";
			}
			else{
				$filename = $callDir.$whichlist[$k];  
			}
			copen(*IN1,"<$filename") or print_usage "Can't open filename $filename. $!";
			print "Reading in the genotyping data from file $filename\n";
			#CATHEDIT: Match column headers with names "SNP Name", "Allele1 - Top", "Allele2 - Top"
			#to work around problem of different geno file batches having different column assignments
			$count = 0;
			while($line=readline(IN1)) {
#				$headarray[$count] = $line;
				if(substr($line, 0, 4) eq "[Dat") { #Finding the line "[Data]" in genotyping file
#
#					$headerline = $count + 1; #Next line contains column headers
					$headerline = readline(IN1); #Next line contains column headers
					last;
				}
#				$count++;
#				if($count>20){		#These 9 hashes dotted about here are to correct a bug.  If all well then can remove dead lines in future.
#					last; 
#				}
			}
#			@header = split(m{\t|\n|\r},$headarray[$headerline]);
			@header = split(m{\t|\n|\r},$headerline);
#			print "HEADER LINE\n".$headarray[$headerline];
			print "HEADER LINE\n$headerline";
			%hash_header_index = array_index(@header);
			$SNP_col = $hash_header_index{"SNP Name"};
			$alleleA_col = $hash_header_index{"Allele1 - Top"};
			$alleleB_col = $hash_header_index{"Allele2 - Top"};
			if(defined($SNP_col) && defined($alleleA_col)  && defined($alleleB_col) ) {
				print "Matched header columns in genotyping file\n";
				print "SNP_col = ".$SNP_col.", alleleA_col = ".$alleleA_col.", alleleB_col = ".$alleleB_col."\n";
			}
			else{
				$SNP_col = 1;
				$alleleA_col = 3;
				$alleleB_col = 4;
				print "Could not match header columns in genotyping file.  Revert to default column indices\n";
			}
			while($line=readline(IN1)){
                @temp=split(m{\t|\n|\r}, $line); #CATHEDIT: Switched from \s+ split to \t split
				if($temp[$SNP_col] =~ /^rs\d+$/ && "$temp[$alleleA_col]$temp[$alleleB_col]" !~ /[ID]/){ # only selecting markers with rsxxx and not will alleles I or D (insertion/deletion)
					if($line_cnt==0){
						print "@temp\n";
					}
					#print "temp[0]=$temp[0],temp[1]=$temp[1]\n";
					# recode data to AA versus A and A etc
					# missing data is -- afterwards.
					#print "temp[4]=$temp[4],temp[5]=$temp[5]\n";
					# Be careful this sometimes needs toggling between 1/3/4 and 2/5/6
					if($temp[$alleleA_col].$temp[$alleleB_col] eq "--") {  #Counting the missing genotypes
						$nocall_cnt[$k]++;
					}
					if(!($temp[$alleleA_col] eq "-") && $temp[$alleleA_col] eq $temp[$alleleB_col]) { #Count homozygous genotypes
						$hom_cnt[$k]++;
					}
					if(!($temp[$alleleA_col] eq $temp[$alleleB_col])) { #Count heterozygous genotypes
						$het_cnt[$k]++;
					}
					$genos_orig{$temp[$SNP_col]}[$k]=$temp[$alleleA_col].$temp[$alleleB_col]; # genos_orig has for each SNP an entry for each genotype in the order of the samples in the whichSamplesFile
					$line_cnt+=1;
					if($line_cnt%100000==0){
						print "Read in $line_cnt genotyping lines of file $filename\n";
					}
				}
			}
			print "Total number of SNPs in the genotyping file with rs designations = $line_cnt\n";
			push(@total_line_cnt,$line_cnt);
			$num_rs_snps[$k] = $line_cnt;
			close(IN1);
		}
		}
		
		#callFile defined for Illumina data containing all sample data...

		else {

			$filename = $callFile;
			copen(*IN1,"<$filename") or print_usage "Can't open filename $filename. $!";

			print "Reading callFile for Illumina data into an array (this may take a minute or two depending on size of callFile)...\n";

			my @callFileIllumina = <IN1>;
			
			close (IN1);

			print "Finished reading callFile into array.  Now looking at the header rows...\n";

			my $numSNPs;
			my $numberHeaderLines = 0;

			$count = 0;

			foreach (@callFileIllumina)	{
			
				++$numberHeaderLines;

#				print "numberHeaderLines = $numberHeaderLines\n";

				if ($_ =~ /Num SNPs\t(\d+)/) {

					$numSNPs = $1;
					print "According to header of $filename, there are $numSNPs SNPs for each sample in this file.\n"; 
				}

				if(substr($_, 0, 4) eq "[Dat") {

						if (!defined($numSNPs)) {

							exit ("Problem with input call file as [Data] line has come without Num SNPs line.");
						}

						++$numberHeaderLines;
						last;
				}
			}

			$headerline = $callFileIllumina[$numberHeaderLines - 1]; #This line should contain column headers

			@header = split(m{\t|\n|\r},$headerline);

			print "Number of header lines: $numberHeaderLines\n";
			print "Header line with column titles:\n$headerline";

			%hash_header_index = array_index(@header);
			$SNP_col = $hash_header_index{"SNP Name"};
			$alleleA_col = $hash_header_index{"Allele1 - Top"};
			$alleleB_col = $hash_header_index{"Allele2 - Top"};
			$sampleName_col = $hash_header_index{"Sample ID"};

			if ( defined($SNP_col) && defined($alleleA_col)  && defined($alleleB_col) && defined($sampleName_col)) {
				print "Matched header columns in genotyping file\n";
				print "SNP_col = ".$SNP_col.", alleleA_col = ".$alleleA_col.", alleleB_col = ".$alleleB_col.", sampleName_col = ".$sampleName_col."\n";
			}
			else {
#				$SNP_col = 1;
#				$alleleA_col = 3;
#				$alleleB_col = 4;
#				print "Could not match header columns in genotyping file.  Revert to default column indices\n";
				die("Unable to match header columns in the genotyping file.\n");
			}

			for ($k=0; $k<=($#which); $k++) {

				print "which for this sample = $which[$k] \n";

				if ($which[$k] == 0) {

					next;
				}

				$line_cnt = 0;

#			$count = 0;

				for ($i = ($numberHeaderLines + $numSNPs * ($which[$k] - 1)); $i <= ($numberHeaderLines - 1 + $numSNPs * $which[$k]); $i++){

					$line = $callFileIllumina[$i];

        	        @temp=split(m{\t|\n|\r}, $line);

					if($temp[$SNP_col] =~ /^rs\d+$/ && "$temp[$alleleA_col]$temp[$alleleB_col]" !~ /[ID]/) { # only selecting markers with rsxxx and not will alleles I or D (insertion/deletion)

						if ($line_cnt == 0) {

							$sampleName = $temp[$sampleName_col];
							print ""."@temp\n";
						}

						if ($sampleName ne $temp[$sampleName_col]) {

							print "Previous sample name: $temp[$sampleName_col]\n";
							print "Current sample name: $sampleName\n";

							die("It appears that either the call file is not sorted correctly (above sample names are not the same), or that the samples do not have all the same SNPs.");
						}

						if($temp[$alleleA_col].$temp[$alleleB_col] eq "--") {  #Counting the missing genotypes

							$nocall_cnt[$k]++;
						}

						if(!($temp[$alleleA_col] eq "-") && $temp[$alleleA_col] eq $temp[$alleleB_col]) { #Count homozygous genotypes

							$hom_cnt[$k]++;
						}

						if(!($temp[$alleleA_col] eq $temp[$alleleB_col])) { #Count heterozygous genotypes

							$het_cnt[$k]++;
						}

						$genos_orig{$temp[$SNP_col]}[$k]=$temp[$alleleA_col].$temp[$alleleB_col]; # genos_orig has for each SNP an entry for each genotype in the order of the samples in the whichSamplesFile
#						$j+=1;		#appears like pointless code, what's it doing here?
						$line_cnt+=1;

						if($line_cnt%100000==0){

							print "Read in $line_cnt genotyping lines of file $filename\n";
						}
					}
				}

				print "Total number of SNPs (excluding INDELs) with \"rs\" SNP IDs in the genotyping file for sample $sampleName = $line_cnt\n";
				push(@total_line_cnt,$line_cnt);
				$num_rs_snps[$k] = $line_cnt;
			}

			print "Clearing Illumina callFile array from computer memory.";

			@callFileIllumina = ();
		}

		foreach $entry (@total_line_cnt){
			$sum+=$entry;
			print "line_cnts=$entry\n";
		}

		if($sum/$line_cnt != ($#which+1)) {

			print "WARNING: Inconsistency with Illumina genotype data between samples.  Samples appear to have unequal number of SNPs (excluding INDELs).\n";
			print "A common reason for this is the actual inclusion of INDEL SNPs for the following reason.\n";
			print "LINKDATAGEN identifies INDEL SNPs by the nature of their alleles (D or I, rather than the usual A, C, G or T).\n";
			print "Therefore, if an individual has missing genotype data for an actual INDEL SNP (-/-), then this INDEL SNP will not be excluded for that sample.\n";

			if (!defined($removeAIS)) {

				print_usage("Check before proceeding, and consider using -removeAIS to override this issue.");
			}
			
			else {

				print "You have opted to use -removeAIS to remove SNP markers that are absent from any Illumina input file.\n";

				print "Now removing SNPs that were not present in all Illumina input files.\n";

				foreach $key (keys %genos_orig) {

					if (scalar @{$genos_orig{$key}} < scalar (@which)) {

						delete($genos_orig{$key});
					}
				}

				print "\nAfter removing SNPs, there are now ", scalar keys %genos_orig, " markers remaining from the Illumina genotype files.\n";
			}
		}
	}

	if($data eq "a" || $data eq "m"){
		$line_cnt=0;
		copen(*IN1,"<$callFile") or print_usage "Can't open genotype call file $callFile. $!"; 
		print "Reading in the genotyping data from file $callFile ...\n";
		if($data eq "a"){
			while($line=readline(IN1)){
				if($line=~/$str/){ # only selecting markers with SNP_A
					@temp=split(/\s+/,$line);
#					$j=0;		#appears like pointless code, why not use $i since it's in the loop??
					for($i=0;$i<=$#which;$i++){
						if(defined($crlmm)){ # now recoding all data to brlmm type
#							$genos_orig{$temp[0]}[$j]=$temp[$which[$i]]-1;		#appears like pointless code, why not use $i since it's in the loop??
							$genos_orig{$temp[0]}[$i]=$temp[$which[$i]]-1;
						}
						else{
							$genos_orig{$temp[0]}[$i]=$temp[$which[$i]];
						}
						#genos_orig has for each SNP an entry for each genotype in the order of the samples in the whichSamplesFile
						if($genos_orig{$temp[0]}[$i] == 0 || $genos_orig{$temp[0]}[$i] == 2) {
							$hom_cnt[$i]++;
						}
						if($genos_orig{$temp[0]}[$i] == 1) {
							$het_cnt[$i]++;
						}
						if($genos_orig{$temp[0]}[$i] == -1) {
							$nocall_cnt[$i]++;
						}
#						$j+=1;		#appears like pointless code, why not use $i since it's in the loop??
					}
					$line_cnt+=1;
					if($line_cnt%100000==0){
						print "Read in $line_cnt genotyping lines\n";
					}
				}
			}
		}
		if($data eq "m"){
			while($line=readline(IN1)){
				@temp=split(/\s+/,$line);
#				$j=0;		#appears like pointless code, why not use $i since it's in the loop??
				for($i=0;$i<=$#which;$i++){
					if(defined($crlmm)){ # now recoding all data to brlmm type
						$genos_orig{$temp[0]}[$i]=$temp[$which[$i]]-1;
					}
					else{
						$genos_orig{$temp[0]}[$i]=$temp[$which[$i]];
					}
					#genos_orig has for each SNP an entry for each genotype in the order of the samples in the whichSamplesFile
					if($genos_orig{$temp[0]}[$i] == 0 || $genos_orig{$temp[0]}[$i] == 2) {
						$hom_cnt[$i]++;
					}
					if($genos_orig{$temp[0]}[$i] == 1) {
						$het_cnt[$i]++;
					}
					if($genos_orig{$temp[0]}[$i] == -1) {
						$nocall_cnt[$i]++;
					}
					
#					$j+=1;		#appears like pointless code, why not use $i since it's in the loop??
				}
				$line_cnt+=1;
				if($line_cnt%100000==0){
					print "Read in $line_cnt genotyping lines\n";
				}
			}
		}
		print "Total number of SNPs in (b/c)rlmm file =$line_cnt\n";
	}

	print "Copying genotype data to memory.\n";

	%genosOrigArchive = %genos_orig;

	open(OUT, ">$tablesDir/genotypeCallRateHomHetMiss.txt") or print_usage("Cannot write file $tablesDir/genotypeCallRateHomHetMiss.txt, $!");

	print "\nGenotype call rates:\n\n";

	print "FID\tIID\tCallRate\tHom\tHet\tMiss\n";
	print OUT "FID\tIID\tCalRate\tHom\tHet\tMiss\n";

	print "===================================\n";

	for($k=0;$k<=($#which);$k++){
		if($data eq "a" || $data eq "m"){
			$den=$line_cnt;
		}
		if($data eq "i"){
			$den=$num_rs_snps[$k];
		}

		print "$whichnames[2*$k]\t$whichnames[2*$k+1]\t";
		print OUT "$whichnames[2*$k]\t$whichnames[2*$k+1]\t";

		$missrate = $nocall_cnt[$k]/$den;
		$homrate = $hom_cnt[$k]/$den;
		$hetrate = $het_cnt[$k]/$den;

		printf OUT ("%1.4f\t%1.4f\t%1.4f\t%1.4f\n", (1-$missrate), $homrate, $hetrate, $missrate);

		printf "%1.4f", (1-$missrate);
		print "\t";
		printf "%1.4f", $homrate;
		print "\t";
		printf "%1.4f", $hetrate;
		print "\t";
		printf "%1.4f", $missrate;
		print "\n";
	}

	close(OUT);

	print "===================================\n";		
}

sub read_in_ped{
my $i;
my $j;

	open(OUT, ">", "$tablesDir/pedigree.txt") or print_usage("Cannot write file $tablesDir/pedigree.txt, $!");

	if($arr3[0]!~/\d+/){
		print "\nWarning: Deleted first line of pedigree file as it is a header line.\n";
		shift(@arr3);
	}

	print "\nFID\tIID\tPID\tMID\tSex\tAffStat\tGenotyped?\n";
	print OUT "FID\tIID\tPID\tMID\tSex\tAffStat\tGenotyped?\n";

	foreach $line (@arr3){
		@temp=split(/\s+/,$line);
		if(@temp != 6){
			print_usage("Error in pedigree file: A line has ".scalar(@temp)." entries instead of 6."); 
		}

		if ( $temp[4] !~ /^[012]{1}$/ ) {

			print_usage("The 5th column of your pedfile should describe the sex of the individuals.\n\nIt can only take values 0 (unknown), 1 (male) or 2 (female). Invalid value identified: $temp[4]"); 		#check sex column for valid values 0, 1 or 2 only
		}

		for($j=0;$j<=$#temp;$j++){
			$ped[$pedno][$j]=$temp[$j];
		}
		$pedno+=1;
	}
	
	if($pedno!=($#whichsamples+1)){
		print_usage("The number of lines in you whichsamples file and your pedigree file differ: whichsamples file = ".($#whichsamples+1)." pedigree file = $pedno");
	}

	for($i=0;$i<$pedno;$i++){

		for($j=0;$j<6;$j++){
			printf ("%.7s\t", $ped[$i][$j]);
			print OUT "$ped[$i][$j]\t";
		}
		if($whichsamples[$i]>0){
			print "Yes\n";
			print OUT "Yes\n";
		}
		else{
			print "No\n";
			print OUT "No\n";
		}
	}

	for($i=0;$i<$pedno;$i++){	#check pedfile for proper construction; parents sex correct and no duplicates in pedfile

		my $motherFound = 0;
		my $fatherFound = 0;

		if (($ped[$i][2] eq 0 && $ped[$i][3] ne 0) || ($ped[$i][2] ne 0 && $ped[$i][3] eq 0)) {

			print_usage("Individuals must have both or neither (0 0) parents defined: $ped[$i][0] $ped[$i][1]'s father = $ped[$i][2] and mother = $ped[$i][3]");
		}

		if ( ( $ped[$i][2] ne 0 && $ped[$i][3] ne 0 ) && ( $ped[$i][2] eq $ped[$i][3]) ) {

			print_usage("Individual's parents must be different: $ped[$i][0] $ped[$i][1]'s father = $ped[$i][2] and mother = $ped[$i][3]");
		}

		if ( ( $ped[$i][2] eq $ped[$i][1] ) ) {

			print_usage("Paternal ID (PID) must be different to individual ID (IID): FID = $ped[$i][0] IID = $ped[$i][1] PID = $ped[$i][2]");
		}

		if ( ( $ped[$i][3] eq $ped[$i][1] ) ) {

			print_usage("Maternal ID (MID) must be different to individual ID (IID): FID = $ped[$i][0] IID = $ped[$i][1] MID = $ped[$i][3]");
		}

		for($j=0;$j<$pedno;$j++){

			if ($i == $j) {
				next;	#identical row in pedfile
			}
			
			if ( ( $ped[$i][0] eq $ped[$j][0] ) && ( $ped[$i][1] eq $ped[$j][1] ) ) {

				print_usage("Duplicate individuals have been found in your pedfile: FID = $ped[$i][0] and IID = $ped[$i][1]");
			}
		}

		if ($ped[$i][2] eq 0 && $ped[$i][3] eq 0) {		#no parents described

			next;
		}

		for($j=0;$j<$pedno;$j++){

			if ($i == $j) {
				next;	#identical row in pedfile
			}
			
			if ( $ped[$i][0] ne $ped[$j][0] ) {

				next;	#different family		
			}
			
			if ( $ped[$i][2] eq $ped[$j][1] ) {		#father found

				if ($ped[$j][4] != 1) {
				
					print_usage("Individual \"$ped[$i][0] $ped[$i][1]\" has a father \"$ped[$j][0] $ped[$j][1]\" whose sex must be labelled male (1) in the pedfile: wrong value = $ped[$j][4]");
				}
				else {

					$fatherFound = 1;	#father found with correct sex
					next;
				}
			}

			if ( $ped[$i][3] eq $ped[$j][1] ) {		#mother found

				if ($ped[$j][4] != 2) {
				
					print_usage("Individual \"$ped[$i][0] $ped[$i][1]\" has a mother \"$ped[$j][0] $ped[$j][1]\" whose sex must be labelled female (2) in the pedfile: wrong value = $ped[$j][4]");
				}
				else {

					$motherFound = 1;	#mother found with correct sex
					next;
				}
			}
		}

		if ($fatherFound == 0) {

			print_usage("Individual \"$ped[$i][0] $ped[$i][1]\" is described as having a father \"$ped[$i][0] $ped[$i][2]\" who can not be found in the pedigree file.\n");
		}

		if ($motherFound == 0) {

			print_usage("Individual \"$ped[$i][0] $ped[$i][1]\" is described as having a mother \"$ped[$i][0] $ped[$i][3]\" who can not be found in the pedigree file.\n");
		}
	}

	close(OUT);
}

#Read in the annotFile, the legacy one (introduced June 2011) with 18 fields:
#Probe_set_ID	rs_name	Chrom	Strand	deCODE_genetic_map_position	physical_position_build37	SNP_type	
#allele_frequencies_CEU	ASW	CHB	CHD	GIH	JPT	LWK	MEX	MKK	TSI	YRI
sub read_in_annot_illumina_legacy{
	my $i;
	my $line;
	my $line_cnt=0;
	my $key;
	my $test=0;
	my $totalline_cnt=0;

	open(OUT3, ">$tablesDir/chrX_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrX_SNPs.txt, $!");
	open(OUT4, ">$tablesDir/chrY_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrY_SNPs.txt, $!");;
	open(OUT5, ">$tablesDir/chrMT_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrMT_SNPs.txt, $!");;

	$line=readline(IN2);
	LINE:while ($line=readline(IN2)){
		$totalline_cnt+=1;
		@temp=split(/\s+/,$line);
		#if($line_cnt < 100) {  #ERROR CHECKING
		#	print STDERR $temp[$pop]."\n";
		#	}
		if(!defined($temp[1])){
			next LINE;
		}
		if($temp[1] =~/[rs]\d+/ && $temp[$popCol] =~ /[\d\.]+/){  #skip SNP if there is no frequency annotation

			if (%regions)           {       #test if regions are being used and if so whether region included...

				my $inRegion = 0;

				if (!defined($regions{$temp[2]})){

					next LINE;
				}

				else {

					foreach my $tempRegion (@{$regions{$temp[2]}}) {

						if ($temp[5] >= $tempRegion->[0] && $temp[5] <= $tempRegion->[1]) {

							$inRegion = 1;

							last;
						}
					}

					if ($inRegion == 0) {

						next LINE;
					}
				}
			}

			if ($temp[6] =~ /[ID]/) {			#This is to exclude 'SNPs' that are defined by insertion/deletion alleles

				next LINE;
			}

			if (!defined($kCT) && $temp[6] =~ /A\/T|T\/A|G\/C|C\/G/) {			#$kCT is an option to over-ride default behaviour here which is to exclude 'SNPs' that are either A/T, T/A, C/G or G/C

				next LINE;
			}

			if (defined($popHetTest) && ($temp[$popCol] < 0.4 || $temp[$popCol] > 0.6)) {	#when running the -popHetTest, this ensures SNPs are chosen within a narrow window of MAF > 0.4

				next LINE;
			}

			if (defined($minMAF) && ($temp[$popCol] < $minMAF || $temp[$popCol] > (1 - $minMAF))) {	#used to specify the minimum MAF acceptable

				next LINE;
			}

			if (defined($maxMAF) && ($temp[$popCol] > $maxMAF && $temp[$popCol] < (1 - $maxMAF))) {	#used to specify the maximum MAF acceptable

				next LINE;
			}

			$annot_orig{$temp[1]}[0]=$temp[0]; # snp name
			$annot_orig{$temp[1]}[1]=$temp[2]; # chromosome
			$annot_orig{$temp[1]}[2]=$temp[6]; # type of SnP eg T/C  **This used to be read in as "[T/C]", now it is "T/C"
			$annot_orig{$temp[1]}[3]=$temp[4]; # sex averaged genetic position
			$annot_orig{$temp[1]}[4]=$temp[$popCol]; # allele frequency of A
			$annot_orig{$temp[1]}[5]=2*$temp[$popCol]*(1-$temp[$popCol]); # heterozygosity
			$annot_orig{$temp[1]}[6]=$temp[0]; # rs name
			$annot_orig{$temp[1]}[7]=$temp[5]; # physical position
			$annot_orig{$temp[1]}[8]=$temp[3]; # strand
			#$annot_orig{$temp[1]}[8]=$temp[40];  #used to be used for TOP and BOT but these field are no longer necessary.
			#$annot_orig{$temp[1]}[9]=$temp[41];
			$line_cnt+=1;	
			if($annot_orig{$temp[0]}[1]=~/X\z/){
				for($i=0;$i<8;$i++){
					print OUT3 "$annot_orig{$temp[0]}[$i]\t";
				}
				print OUT3 "\n";
			}
			if($annot_orig{$temp[0]}[1]=~/Y\z/){
				for($i=0;$i<8;$i++){
					print OUT4 "$annot_orig{$temp[0]}[$i]\t";
				}
				print OUT4 "\n";
			}
			if($annot_orig{$temp[0]}[1]=~/M\z/){
				for($i=0;$i<8;$i++){
					print OUT5 "$annot_orig{$temp[0]}[$i]\t";
				}
				print OUT5 "\n";
			}
			if($line_cnt%100000==0){
				print "$line_cnt SNPs so far selected from the annotation file...\n";
			}
		}
	}	
		
	close(OUT3);
	close(OUT4);
	close(OUT5);
	print "Total number of SNPs in the Illumina map annotation file = $totalline_cnt\n";
#	print "Total number of usable SNPs in Illumina annotation file from column $pops[$popCol - $popColOriginal] = $line_cnt\n";

if (defined($bestPopTest)) {
	print "Total number of usable SNPs in Illumina annotation file from column $pops[$popCol - $popColOriginal] = $line_cnt\n";
} elsif (defined($pop)) {
	print "Total number of usable SNPs in Illumina annotation file from column $pop = $line_cnt\n";
#} elsif (defined($pop)) {
#	print "Total number of usable SNPs in Illumina annotation file from column " . ($pop eq "USR" ? $popCol + 1 : $pop ) . " = $line_cnt\n";
} else {
	print "Total number of usable SNPs in Illumina annotation file from CEU column (default selection) = $line_cnt\n";
}

}

#Chrom	
#physical_position_(bp)	
#rs_name	
#Strand	
#TOPBOT strand
#deCODE_genetic_map_position	
#alleleA/alleleB
#allele_frequencies_1)CEU 2)ASW	3)CHB 4)CHD 5) GIH 6)JPT 7)LWK 8)MEX 9)MKK10)TSI 10)YRI

# Chrom physical_position_build37 rs_name Strand TOPBOT deCODE_genetic_map_position A B allele_frequencies_CEU ASW CHB CHD GIH JPT LWK MEX MKK TSI YRI
sub read_in_annot_mps_illuminaSNP{
my $i;
my $line;
my $key;
my $test=0;
my $totalline_cnt=0;
my $line_cnt=-1;
my $head_comments = 1;
my $annot_version_check = 1;
my $rsname_col; 
my $chr_col;
my $strand_col;
my $topbot_col;
my $genetic_pos_col;
my $physical_pos_col;
my $freq_col;
my $alleleA_col;
my $alleleB_col;
my @header_array;

	open(OUT3, ">$tablesDir/chrX_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrX_SNPs.txt, $!");
	open(OUT4, ">$tablesDir/chrY_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrY_SNPs.txt, $!");;
	open(OUT5, ">$tablesDir/chrMT_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrMT_SNPs.txt, $!");;

	LINE:while ($line=readline(IN2)){
        chomp $line;
		$totalline_cnt+=1;
        if($head_comments) {
            # Check head comment lines
            if($line =~ /^#/) {
                # For old version format
                if($line =~ /^# Minimum LINKDATAGEN revision: (\d+)$/) {
                    if(997 < $1) {
                        # this script is too old
                        print_usage("Annotation file is requires newer version of LINKDATAGEN. Make sure you have not changed the annotation header or contact the developers. (SVN Revision $1 required, but have $script_revision)"); 
                    }
                    $annot_version_check = 1;
                }
                # For new version format:
                my @min_ldg_for_annot;
                if(@min_ldg_for_annot = $line =~ /^# Minimum LINKDATAGEN revision: v(\d+)\.(\d+)\.(\d+)$/) {
                    $script_revision =~ /^(\d+)\.(\d+)\.(\d+)$/;
                    if($1 < $min_ldg_for_annot[0] || $2 < $min_ldg_for_annot[1] || $3 < $min_ldg_for_annot[2]) {
                        # this script is too old
                        $line =~ /^# Minimum LINKDATAGEN revision: v(\d+\.\d+\.\d+)$/;
                        print_usage("Annotation file is requires newer revision of LINKDATAGEN. (Version $1 required, but have Perl script v$script_revision)"); 
                    }
                    $annot_version_check = 1;
                }
                next LINE;
            } else {
                $head_comments = 0;
                unless($annot_version_check) {
                    print_usage("Annotation minimum version not found. Are you using the latest annotation file?");
                }
            }
        }

		@temp=split(/\s+/,$line);

        if($line_cnt == -1) {
            @header_array = @temp;
            my %hash_header_index = array_index(@header_array);
            $rsname_col = $hash_header_index{"rs_name"};
            $chr_col = $hash_header_index{"Chrom"};
            $strand_col = $hash_header_index{"Strand"};
            $topbot_col = $hash_header_index{"TOPBOT"};
            $genetic_pos_col = $hash_header_index{"deCODE_genetic_map_position"};
            $physical_pos_col = $hash_header_index{"physical_position_build37"} || $hash_header_index{"physical_position_build38"};
            if($pop eq "USR" || defined($bestPopTest)) {
                # User defined column
                $freq_col = $popCol;
            } else {
                $freq_col = defined($hash_header_index{"$pop"}) ? $hash_header_index{"$pop"} : $hash_header_index{"allele_frequencies_$pop"}; 
                if($popCol != $freq_col) {
                    die "Internal error (contact developers), population column isn't as expected";
                }
            } 
            $alleleA_col = $hash_header_index{"A"};
            $alleleB_col = $hash_header_index{"B"};
            if(defined($rsname_col) &&  defined($chr_col) &&  defined($strand_col) &&  defined($genetic_pos_col) &&  defined($physical_pos_col) &&  
                defined($freq_col) && defined($topbot_col) && defined($alleleA_col) && defined($alleleB_col) ){
                print "Matched column indices in annotation file\n";
            }
            else{
                print_usage("Unable to match column indices in annotation file. The file may be corrupted.\n" . 
                    "SNP_name => $rsname_col, Chromosome => $chr_col, Strand => $strand_col, TOP/BOT => $topbot_col, Genetic_position => $genetic_pos_col, Physical_position => $physical_pos_col, alleleA_col => $alleleA_col, alleleB_col => $alleleB_col, Frequency => $freq_col");
            }
            $line_cnt++;
            next LINE;
        }

        $temp[$chr_col] =~ s/^chr//; # Allow the newer annotation to continue working properly without a massive recode

		if($temp[$freq_col] =~/\d+/){

			if (%regions)           {       #test if regions are being used and if so whether region included...
				my $inRegion = 0;
				if (!defined($regions{$temp[$chr_col]})){
					next LINE;
				}

				else {

					foreach my $tempRegion (@{$regions{$temp[$chr_col]}}) {
						if ($temp[$physical_pos_col] >= $tempRegion->[0] && $temp[$physical_pos_col] <= $tempRegion->[1]) {
							$inRegion = 1;
							last;
						}
					}

					if ($inRegion == 0) {
                        # We are not in region for output
						next LINE;
					}
				}
			}

			if (defined($popHetTest) && ($temp[$freq_col] < 0.4 || $temp[$freq_col] > 0.6)) {	#when running the -popHetTest, this ensures SNPs are chosen within a narrow window of MAF > 0.4
				next LINE;
			}

			if (defined($minMAF) && ($temp[$freq_col] < $minMAF || $temp[$freq_col] > (1 - $minMAF))) {	#used to specify the minimum MAF acceptable
				next LINE;
			}

			if (defined($maxMAF) && ($temp[$freq_col] > $maxMAF && $temp[$freq_col] < (1 - $maxMAF))) {	#used to specify the maximum MAF acceptable
				next LINE;
			}

            my $snp_name = $temp[$rsname_col];

            if ($temp[$strand_col] eq '-') {
                print STDERR "Warning: annotation file contains SNPs on the reverse strand. Not readable at this time. Line:\t";
                for($i=0;$i<ANNOT_TOTAL_COLS;$i++){
					print STDERR "$annot_orig{$snp_name}[$i]\t";
				}
				print STDERR "\n";
				next LINE;
            }

			$annot_orig{$snp_name}[ANNOT_SNP] = $snp_name; # snp name
			$annot_orig{$snp_name}[ANNOT_CHR] = $temp[$chr_col];  # chromosome
			$annot_orig{$snp_name}[ANNOT_ALA] = $temp[$alleleA_col]; # allele A	
			$annot_orig{$snp_name}[ANNOT_ALB] = $temp[$alleleB_col]; # allele B	
			$annot_orig{$snp_name}[ANNOT_GPO] = $temp[$genetic_pos_col]; # sex averaged genetic position
			$annot_orig{$snp_name}[ANNOT_FRE] = $temp[$freq_col]; # allele frequency of A
			$annot_orig{$snp_name}[ANNOT_HET] = heterozygosity_from_af($temp[$freq_col]); #heterozygosity
			$annot_orig{$snp_name}[ANNOT_RSN] = $snp_name; # rs name
			$annot_orig{$snp_name}[ANNOT_POS] = $temp[$physical_pos_col]; # physical position
			$annot_orig{$snp_name}[ANNOT_TBO] = $temp[$topbot_col]; # Reference TOP/BOT strand

			$line_cnt+=1;
			if($annot_orig{$snp_name}[ANNOT_CHR]=~/X\z/){
				for($i=0;$i<ANNOT_TOTAL_COLS;$i++){
					print OUT3 "$annot_orig{$snp_name}[$i]\t";
				}
				print OUT3 "\n";
			}
			if($annot_orig{$snp_name}[ANNOT_CHR]=~/Y\z/){
				for($i=0;$i<ANNOT_TOTAL_COLS;$i++){
					print OUT4 "$annot_orig{$snp_name}[$i]\t";
				}
				print OUT4 "\n";
			}
			if($annot_orig{$snp_name}[ANNOT_CHR]=~/M\z/){
				for($i=0;$i<ANNOT_TOTAL_COLS;$i++){
					print OUT5 "$annot_orig{$snp_name}[$i]\t";
				}
				print OUT5 "\n";
			}
			if($line_cnt%100000==0){
				print "$line_cnt SNPs so far selected from the annotation file...\n";
			}
		}
	}

	close(OUT3);
	close(OUT4);
	close(OUT5);
	close(IN2);

	print "Total number of SNPs in the MPS Hapmap annotation file = $totalline_cnt\n";
	print "Total number of usable SNPs in the MPS Hapmap annotation file from column $header_array[$freq_col] = $line_cnt\n";
	
	return ($line_cnt);
}

# read in annotation files are now taken from our internal format as generated by Catherine.
#Affy_SNP_name	
#rs_name	
#Chrom	
#Strand	
#deCODE_genetic_map_position	
#physical_position_(bp)	
#allele_frequencies_1)CEU 2)ASW	3)CHB 4)CHD 5) GIH 6)JPT 7)LWK 8)MEX 9)MKK10)TSI 10)YRI
sub read_in_annot_affy{
my $i;
my $line;
my $line_cnt=0;
my $key;
my $test=0;
my $totalline_cnt=0;
	# get rid of header line	
	$line=readline(IN2);

	open(OUT3, ">$tablesDir/chrX_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrX_SNPs.txt, $!");;
	open(OUT4, ">$tablesDir/chrY_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrY_SNPs.txt, $!");;
	open(OUT5, ">$tablesDir/chrMT_SNPs.txt") or print_usage("Cannot write file $tablesDir/chrMT_SNPs.txt, $!");;

	LINE:while ($line=readline(IN2)){
		$totalline_cnt+=1;
		@temp=split(/\s+/,$line);
		if($temp[$popCol] =~/\d+/ && $temp[2] =~/\w+/ && $temp[4] =~/\w+/ && $temp[5] =~/\w+/){
			#June 2011: extra conditions on the if statement.  The SNP will not be read into annot_orig if
			#it has missing map info (chr, genetic map or physical map)

			if (%regions)           {       #test if regions are being used and if so whether region included...

				my $inRegion = 0;

				if (!defined($regions{$temp[2]})){

					next LINE;
				}

				else {

					foreach my $tempRegion (@{$regions{$temp[2]}}) {

						if ($temp[5] >= $tempRegion->[0] && $temp[5] <= $tempRegion->[1]) {

							$inRegion = 1;

							last;
						}
					}

					if ($inRegion == 0) {

						next LINE;
					}
				}
			}

			if (defined($popHetTest) && ($temp[$popCol] < 0.4 || $temp[$popCol] > 0.6)) {	#when running the -popHetTest, this ensures SNPs are chosen within a narrow window of MAF > 0.4

				next LINE;
			}

			if (defined($minMAF) && ($temp[$popCol] < $minMAF || $temp[$popCol] > (1 - $minMAF))) {	#used to specify the minimum MAF acceptable

				next LINE;
			}

			if (defined($maxMAF) && ($temp[$popCol] > $maxMAF && $temp[$popCol] < (1 - $maxMAF))) {	#used to specify the maximum MAF acceptable

				next LINE;
			}

			$annot_orig{$temp[0]}[0]=$temp[0]; # snp name
			$annot_orig{$temp[0]}[1]=$temp[2]; # chromosome
			$annot_orig{$temp[0]}[2]=$temp[9]; # allele A
			$annot_orig{$temp[0]}[3]=$temp[4]; # sex averaged genetic position
			$annot_orig{$temp[0]}[4]=$temp[$popCol]; # allele frequency of A
			$annot_orig{$temp[0]}[5]=2*$temp[$popCol]*(1-$temp[$popCol]); # heterozygosity
			$annot_orig{$temp[0]}[6]=$temp[1]; # rs name
			$annot_orig{$temp[0]}[7]=$temp[5]; # physical position
			$line_cnt+=1;	
			if($annot_orig{$temp[0]}[1]=~/X\z/){
				for($i=0;$i<8;$i++){
					print OUT3 "$annot_orig{$temp[0]}[$i]\t";
				}
				print OUT3 "\n";
			}
			if($annot_orig{$temp[0]}[1]=~/Y\z/){
				for($i=0;$i<8;$i++){
					print OUT4 "$annot_orig{$temp[0]}[$i]\t";
				}
				print OUT4 "\n";
			}
			if($annot_orig{$temp[0]}[1]=~/M\z/){
				for($i=0;$i<8;$i++){
					print OUT5 "$annot_orig{$temp[0]}[$i]\t";
				}
				print OUT5 "\n";
			}
			if($line_cnt%100000==0){
				print "$line_cnt SNPs so far selected from the annotation file...\n";
			}
		}
	}
	close(OUT3);
	close(OUT4);
	close(OUT5);
	print "Total number of SNPs in the Affy annotation file = $totalline_cnt\n";
	print "Total number of usable SNPs in Affymetrix annotation file from column $pops[$popCol - $popColOriginal] = $line_cnt\n";
}

sub removeSNPs {

	my $line;

	while($line=readline(FRS)){

		chomp($line);
		delete($annot_orig{$line});
	}
	
	print "\nAfter removing SNPs, there are now ", scalar keys %annot_orig, " markers remaining from the annotation file.\n";
}

sub keepSNPs {

	my $line;
	my %keepSNPs;

	while($line=readline(FKS)){

		chomp($line);
		$keepSNPs{$line} = 1;
	}

	foreach $line (keys %annot_orig){

		if(!defined($keepSNPs{$line})){

			delete($annot_orig{$line});
		}
	}
	
	print "\nAfter keeping SNPs, there are now ", scalar keys %annot_orig, " markers remaining from the annotation file.\n";
}

# check that all markers are covered by both the genotyping data and the annotation files.
sub check_geno_markers{

my $key;
my $cnt=0;
my @to_remove=();
	print "\nChecking the markers against genotype data and annotation file data and vice versa.\n";
	print "Markers absent in either cannot be included in the final datasets.\n";

	print "\nThere are currently:\n\n";
	print "\t", scalar keys %annot_orig, " markers from the annotation file.\n";
	print "\t", scalar keys %genos_orig, " markers from the genotype file(s).\n";
	
	foreach $key (keys %genos_orig){
		if(!defined($annot_orig{$key})){
			#print "SNP $key is in genotyping data file and not in the annotation file\n";
			$cnt+=1;
			# remove this SNP from the genotyping data file.
			push(@to_remove,$key);
		}
	}
	
	print "\nTotal number of markers present in the genotype file(s) but absent from the annotation file: ", scalar @to_remove;
	print "\nRemoving these markers from the genotyping data...\n";

	foreach $key (@to_remove){
		delete($genos_orig{$key});
	}

	# check no of lines now:

	@to_remove=();

	foreach $key (keys %annot_orig){
		if(!defined($genos_orig{$key})){
			push(@to_remove,$key);
		}
	}
	
	print "\nTotal number of markers present in the annotation file but absent from the genotype file(s): ", scalar @to_remove;
	print "\nRemoving these markers from the annotation data...\n";

	foreach $key (@to_remove){
		delete($annot_orig{$key});
	}

	print "\nThere are currently:\n\n";
	print "\t", scalar keys %annot_orig, " markers from the annotation file.\n";
	print "\t", scalar keys %genos_orig, " markers from the genotype file(s).\n";
}

# remove any "within-family homozygosity-by-state" SNP markers.  can be intersection or union across families.
sub removeWFHBS{
my $key;
my $i;
my $m;
my $test;
my %reject_snps=();
my %reject_snps_byfamily=();
my $k;
my $keycnt=0;
my $testvalue=2;
#my %fam=();
#my %fam_i=();
my %famgeno=();
#my %famgenocnt=();
#my $family;
#my $entry;
my $test1;
my $test2;
#my $familiesTotal;
my $removed = 0;

## Much of this SORTING code is now redundant as has been placed higher up in the code, to stop duplicity.
#	#sort  data into families
#	for($i=0;$i<=$#whichsamples;$i++){
#		$fam{$ped[$i][0]}{$ped[$i][1]}=$i;
#	}

#	$familiesTotal = scalar keys (%fam);

	foreach $family (keys %fam){
		$reject_snps_byfamily{$family}=0;
	}
		#print "family $family\n";
#		foreach $entry (keys %{$fam{$family}}){
#			if($whichsamples[$fam{$family}{$entry}]!=0){
				#print "Individual $entry is pedentry $fam{$family}{$entry}, with which=$whichsamples[$fam{$family}{$entry}]\n";
#				for($i=0;$i<=$#which;$i++){
#					if($whichsamples[$fam{$family}{$entry}]==$which[$i]){
#						push(@{$fam_i{$family}},$i);
#						last;
#					}				
#				}	#push(@{$fam_i{$family}},$whichsamples[$fam{$family}{$entry}]);
#				$famgenocnt{$family}+=1;
#			}
#		}
#	}
	# need to perform WFHBS testing now w/n each family separately
	# limit which $i to run through. Easiest done by changing fam to hash/hash with $i as value?
	# use 
	SNP:foreach $key (keys %genos_orig){ 
		foreach $family (keys %fam){
		
			if (defined($famgenocnt{$family})) {
#			if($keycnt<=$testvalue){
#				print "Family $family\n";
#			}
			$famgeno{$family}{"-1"}=0;
			$famgeno{$family}{"0"}=0;
			$famgeno{$family}{"1"}=0;
			$famgeno{$family}{"2"}=0;
			for($i=0;$i<=$#{$fam_i{$family}};$i++){
				$famgeno{$family}{$genos_orig{$key}[$fam_i{$family}[$i]]}+=1;
#				if($keycnt<=$testvalue){
#					print "genos_orig{$key}[$fam_i{$family}[$i]]=$genos_orig{$key}[$fam_i{$family}[$i]]\n";
#				}
			}
			# tests for WFHBS
			$test1=$famgeno{$family}{"-1"}+$famgeno{$family}{"0"};
			$test2=$famgeno{$family}{"-1"}+$famgeno{$family}{"2"};
			if($test1==$famgenocnt{$family} ||$test2==$famgenocnt{$family}){
				if($keycnt<=$testvalue){
					print "$key is a \"within-family homozygosity-by-state\" SNP marker in family $family\n";
				}
				$reject_snps{$key}{$family}=1;
				$reject_snps_byfamily{$family}+=1;
			}
		}
		else{
				$reject_snps{$key}{$family}=1;
				$reject_snps_byfamily{$family}+=1;
		}
		}
		$keycnt+=1;
	}
	print "\nNumber of snps prior to \"within-family homozygosity-by-state\" SNP marker removal = ",scalar keys %annot_orig,"\n\n";
	print "Family\t#WFHBS SNP Markers\n";
	print "============================================\n";
	foreach $family (keys %fam){
		print "$family\t",$reject_snps_byfamily{$family},"\n";
	}

	print "\nTotal of number of families = $familiesTotal\n";

	print "\nUnion of \"within-family homozygosity-by-state\" SNP markers across all families = ",scalar keys (%reject_snps),"\n\n";

	if ($removeWFHBS eq "u") {

		print "Removing SNP markers that are \"within-family homozygosity-by-state\" in ANY family - that is the union.\n";

		foreach $key (keys %reject_snps){

			delete($genos_orig{$key});
			delete($annot_orig{$key});

			$removed++;
		}

		print "\nNumber of \"within-family homozygosity-by-state\" SNP markers removed (from union across families) = $removed\n";
	}

	elsif ($removeWFHBS eq "i") {

		print "Removing SNP markers that are \"within-family homozygosity-by-state\" in ALL families - that is the intersection.\n";

		foreach $key (keys (%reject_snps)) {

			if (scalar keys (%{$reject_snps{$key}}) == $familiesTotal) {

				delete($genos_orig{$key});
				delete($annot_orig{$key});

				$removed++;
			}
		}

		print "\nNumber of \"within-family homozygosity-by-state\" SNP markers removed (from intersection across families) = $removed\n";
	}

	else {

		print_usage("Invalid option for -removeWFHBS.\n");
	}

	print "\nThere are currently:\n\n";
	print "\t", scalar keys %annot_orig, " markers from the annotation file.\n";
	print "\t", scalar keys %genos_orig, " markers from the genotype file(s).\n";

}

# extra function for illumina data to make genotypes in the AA etc format into the -1, 0, 1, 2 brlmm format
# assume Rust's setup: the given frequency is the frequency of the alphabetically LOWER allele
sub recode_illumina_data{
my $temp;
#my $cnter=0;	## I am sure this is defunct code, and doesn't do anything useful at all anyway. Tom 12/11/12
my $first;
my $second;
my $i;
my $new;
my %snp=(); 
my %strand=();
my $cnt_miss_snps=0;
my $topbot;
my %recgeno=();# double hash for recoding purposes
my %flipper=( 
    A => 'T', 
    T => 'A',
    G => 'C', 
    C => 'G',
);
my $m;

	foreach my $snp_name (keys %annot_orig){
		#print "annot_orig{$snp_name}[8]=$annot_orig{$snp_name}[8]\n";
        if($annotFileIlluminaType eq 'old_illumina') {
            if($annot_orig{$snp_name}[2]=~/([ACTG])\/([ACGT])/){
                $first=$1;
                $second=$2;
            }
            # Assuming unambiguous SNPs only, this is only for the old annotation files
            # A becomes the lower allele
            # eg AA, AT, TT -> 0, 1, 2
            if(ord($first) > ord($second) ) { # do the flip of strands to match
                $first = $flipper{$first};
                $second = $flipper{$second};
            }
            if($first ne 'A') {
                die "Internal error, contact developers. Unambiguous assumption broken for legacy Illumina files for SNP $snp_name";
            }
        } else {
            $first  = $annot_orig{$snp_name}[ANNOT_ALA] or die "Internal error, no allele A in annotation at $snp_name";
            $second = $annot_orig{$snp_name}[ANNOT_ALB] or die "Internal error, no allele B in annotation at $snp_name";
            $topbot = $annot_orig{$snp_name}[ANNOT_TBO] or die "Internal error, no TOP/BOT strand in annotation at $snp_name";
            # If on the TOP strand in the annotation, then no need to change
            # If on the BOT strand, need to reverse annotation
            if($topbot eq 'B') {
                # flip the strand to match
                $first = $flipper{$first};
                $second = $flipper{$second};	
                $annot_orig{$snp_name}[ANNOT_FRE] = 1 - $annot_orig{$snp_name}[ANNOT_FRE];
            } 
        }
        if($first eq $second) {
            die "Internal error, contact developers. Illumina recode has first and second allele the same for $snp_name with allele $first";
        }
		if($first gt $second){
            ($first, $second) = ($second, $first);
        }
        $snp{$snp_name} = $first; # store what is the lower allele.
        $recgeno{$snp_name}{"$first$first"}   = 0;
        $recgeno{$snp_name}{"$first$second"}  = 1;
        $recgeno{$snp_name}{"$second$second"} = 2;
        $recgeno{$snp_name}{"--"}              = -1;
	}

	foreach my $snp_name (keys %genos_orig){
	
		if ($genos_orig{$snp_name}[0] eq -1 || $genos_orig{$snp_name}[0] eq 0 || $genos_orig{$snp_name}[0] eq 1 || $genos_orig{$snp_name}[0] eq 2) {
			next;	#simple test to see if it has been recoded already.  linkdatagen needs re-coding to stop this circle of events.
		}

		for($i=0;$i<=$#which;$i++){
            no warnings; # Checking if defined later
			$temp=$recgeno{$snp_name}{$genos_orig{$snp_name}[$i]};
            use warnings;
			if(defined($temp)){
				$genos_orig{$snp_name}[$i]=$temp;
			}
			else{ # missing genotype data
				$genos_orig{$snp_name}[$i] = -1;
				$cnt_miss_snps+=1;
			}
		}
	}
	print "\nTotal number of genotypes for which we could not recode properly: $cnt_miss_snps\n";
    print "\nNote that these $cnt_miss_snps genotypes will be converted to \"missing\" by LINKDATAGEN.\n";
}

sub print_data_test{
my $j;
my $keycnt=0;
my $key;
	foreach $key (keys %annot_orig){
		print "SNP = $key\n";
		for($j=0;$j<=$#which;$j++){
			print "$genos_orig{$key}[$j] ";		
		}
		$keycnt+=1;
		if($keycnt==10){
			last;
		}
		print "\n";
	}
}

# find all Mendelian error SNPs; keep these SNPs in the annotation and genotype hashes if -keepME is applied - default is to remove SNPs.
sub mendelian_errors{
my $key;
my $i;
my $j;
my $k;
my $m;
my $l;
my $c;
my $chr;
my $chromtottot=0;
my $mytot=0;
my $genom;
my $genof;
my $genoc;
my %reject_snps=();
my $sum_total=0;
my %xchromsnpcnt=();
my %ychromsnpcnt=(); #CATHEDIT
my $yprop;
my %xcheck=();
my %chromsnpcnttot=();
my %chromsnpcnt=();
my $females=0;
my $malemean=0;
my $femalemean=0;
my %sexcheck=();
my %ycheck = (); #CATHEDIT
my %zstat=();
my $chrom=0;
my $ychrom=0;
my @chrom_tot=();
my @chrom_mend=();
my $males=0;
my %ind=();
my $sum_inds=0;
my $thisone;
my %check_poss=();
my %snp_cnts=();
my $totmen=0;
my $mendelian_snp=0;
my @chr_order=(); # array of arrays with SNP names as entries
my $whichsnp;
my $whichcM;
my @snporder_sorted=();

my %FSME=(); ### Family Specific Mendelian Errors

	open(OUTM, ">$tablesDir/MendelErrors.txt") or print_usage "Can't write file $tablesDir/MendelErrors.txt, $!\n";
	foreach $key (keys %genos_orig){ 
		$reject_snps{$key}=0;
	}
	for($i=0;$i<$pedno;$i++){ 
		$chromsnpcnttot{$i}=0;
		for($c=0;$c<=($numchrom-1);$c++){
			$chr=$c+1;
			if($c==22){
				$chr="X";
			}
			if($c==23){
				$chr="Y";
			}
			$chromsnpcnt{$i}{$chr}=0;
		}
	}
	# sex check first....
	for($i=0;$i<$pedno;$i++){ # individual to be tested
		$sexcheck{$ped[$i][0].$ped[$i][1]}=0;
		$ycheck{$ped[$i][0].$ped[$i][1]}=0;
		$xcheck{$ped[$i][0].$ped[$i][1]}=0;
		$xchromsnpcnt{$ped[$i][0].$ped[$i][1]}=0;
		$ychromsnpcnt{$ped[$i][0].$ped[$i][1]}=0;
		$ind{$i}=0;
		if($whichsamples[$i]!=0){ # individual has been genotyped
			for($m=0;$m<=$#which;$m++){
				if($which[$m]==$whichsamples[$i]){
					$thisone=$m;	#m is nth (staring from 0) genotyped sample from the whichsamples file... i.e. for 0 0 0 12 0 3  then m=0 is for 12, m=1 for 3 etc...
					last;
				}
			}
			#print "ped[$i][1]=$ped[$i][0].$ped[$i][1],whichsamples[$i]=$whichsamples[$i],which[$m]=$thisone\n";
			#should the ind finder be $which[$m] instead?\n";
			foreach $key (keys %annot_orig){
				if(!defined($annot_orig{$key}[1])){
					print_usage("annot_orig{$key}[1]=$annot_orig{$key}[1]");
				}
				if($annot_orig{$key}[1]=~/\d+/){
					if($genos_orig{$key}[$thisone] > -1){
						$chromsnpcnt{$i}{$annot_orig{$key}[1]}+=1;
					}
				}
				if($annot_orig{$key}[1] eq "Y"){ # Checking markers on Y chromosome
					$ychromsnpcnt{$ped[$i][0].$ped[$i][1]}+=1; # this def is slightly diff to the X def of the same, count all y snps, will be the same for all inds.
					if($genos_orig{$key}[$thisone] > -1){
						#print "called SNP ".$key." genotype ".$genos_orig{$key}[$thisone]."\n";
						$ycheck{$ped[$i][0].$ped[$i][1]}+=1; # count SNPs with called alleles on Y chr
					}
				}
				if($annot_orig{$key}[1]=~/X\z/){ # X chromosome needs special checks
					if($genos_orig{$key}[$thisone] > -1){
						$xchromsnpcnt{$ped[$i][0].$ped[$i][1]}+=1; # no of genotyped x snps
						$chromsnpcnt{$i}{$annot_orig{$key}[1]}+=1;
					}
					if($genos_orig{$key}[$thisone]==1){
						$sexcheck{$ped[$i][0].$ped[$i][1]}+=1; # counts the no of hets
						if($ped[$i][4]==1){
							$genos_orig{$key}[$thisone]= -1; # fixing males
							$reject_snps{$key}+=1;
							$FSME{$ped[$i][0]}{$key}=1;
						}
					}
				}
			}
			$ind{$i}=$sexcheck{$ped[$i][0].$ped[$i][1]}; # only keep male hets as these are errors.
			#print "\ni=$i, which[$m]=$which[$m], thisone= $thisone, sexcheck{$ped[$i][0].$ped[$i][1]}=$sexcheck{$ped[$i][0].$ped[$i][1]}\n";
			#print "xchromsnpcnt=$xchromsnpcnt{$ped[$i][0].$ped[$i][1]}\n";
			if($ped[$i][4]==2){
				$ind{$i}=0; # only keep male hets as these are errors.
			}
			if($xchromsnpcnt{$ped[$i][0].$ped[$i][1]}>0){ # needed a special case for males with no x 
				$sexcheck{$ped[$i][0].$ped[$i][1]}/=$xchromsnpcnt{$ped[$i][0].$ped[$i][1]};
			}
			else{
				undef($sexcheck{$ped[$i][0].$ped[$i][1]});
			}
		}
	}
	$males=0;
	$females=0;
	for($i=0;$i<$pedno;$i++){ # individual to be tested by looking for his/her parents and checking for mendelian errors
		if($whichsamples[$i]!=0){ # individual has been genotyped
			$ychrom=$ychromsnpcnt{$ped[$i][0].$ped[$i][1]};
			if($ped[$i][4]==1){
				if(defined($sexcheck{$ped[$i][0].$ped[$i][1]})){
					$malemean+=$sexcheck{$ped[$i][0].$ped[$i][1]};
					$males+=1;
				}
			}
			if($ped[$i][4]==2){
				if(defined($sexcheck{$ped[$i][0].$ped[$i][1]})){
					$femalemean+=$sexcheck{$ped[$i][0].$ped[$i][1]};
					$females+=1;
				}
			}
		}
	}
	if($males>0){
		$malemean/=$males;
	}
	if($females>0){
		$femalemean/=$females;
	}
	# normalised N(0,1) proportions
	if($males > 0) {
		for($i=0;$i<$pedno;$i++){ # individual to be tested by looking for his/her parents and checking for mendelian errors		
			if($whichsamples[$i]!=0){ # individual has been genotyped
				if(defined($sexcheck{$ped[$i][0].$ped[$i][1]}) && $ped[$i][4]==1){
					if ($malemean==0){
						$zstat{$ped[$i][0].$ped[$i][1]} = 0;
					}
					else{
						$zstat{$ped[$i][0].$ped[$i][1]}=($sexcheck{$ped[$i][0].$ped[$i][1]}-$malemean)/(sqrt($malemean*(1-$malemean)/$males));  
					}
				}
			}
		}
	}
	if($females > 0) {
		for($i=0;$i<$pedno;$i++){ # individual to be tested by looking for his/her parents and checking for mendelian errors		
			if($whichsamples[$i]!=0){ # individual has been genotyped
				if(defined($sexcheck{$ped[$i][0].$ped[$i][1]}) && $ped[$i][4]==2){
					if ($femalemean==0){
						$zstat{$ped[$i][0].$ped[$i][1]} = 0;
					}
					else {
						$zstat{$ped[$i][0].$ped[$i][1]}=($sexcheck{$ped[$i][0].$ped[$i][1]}-$femalemean)/(sqrt($femalemean*(1-$femalemean)/$females));
					}
				}
			}
		}
	}
	print "\nThere are ",$males," genotyped male(s) and ",$females," genotyped female(s).\n";
	print "\nMendelian Error distribution on the X chromosome based only on observed heterozygotes per genotyped individual:\n\n";

	open(OUT, ">$tablesDir/chrX_HetRatePerSample.txt") or print_usage("Cannot write file $tablesDir/chrX_HetRatePerSample.txt, $!");

	print "Pedigree\tIndividual\tSex\tProportion of hets\tSex Specific Zstat\n";
	print OUT "FID\tIID\tSex\tProHetX\tSex_Specific_z-Stat\n";

	print "===========================\n";

	for($i=0;$i<$pedno;$i++){ # individual to be tested by looking for his/her parents and checking for mendelian errors

		print "$ped[$i][0]\t$ped[$i][1]\t";
		print OUT "$ped[$i][0]\t$ped[$i][1]\t";

		if($ped[$i][4]==1){ # male
			print "male\t";
			print OUT "male\t";
		}

		if($ped[$i][4]==2){ #female
			print "female\t";
			print OUT "female\t";
		}

		if($ped[$i][4]==0){ #unknown
			print "unknown\t";
			print OUT "unknown\t";
		}

		if($whichsamples[$i]!=0){ # individual has been genotyped

			if(defined($sexcheck{$ped[$i][0].$ped[$i][1]})){
				printf "%1.4f",$sexcheck{$ped[$i][0].$ped[$i][1]};
				print "\t";
				printf OUT ("%1.4f\t", $sexcheck{$ped[$i][0].$ped[$i][1]});
			}
			else {

				print "na\t";
				print OUT "na\t";
			}

			if(defined($zstat{$ped[$i][0].$ped[$i][1]})){
				printf "%1.4f",$zstat{$ped[$i][0].$ped[$i][1]};
				printf OUT ("%1.4f", $zstat{$ped[$i][0].$ped[$i][1]});
			}
			else {
				print "na";
				print OUT "na";
			}
		}

		else {

			print "na\tna";
			print OUT "na\tna";
		}

		print "\n";
		print OUT "\n";
	}

	close(OUT);

	print "===========================\n\n";

	open(OUT, ">$tablesDir/chrY_CallRatePerSample.txt") or print_usage "Cannot write file $tablesDir/chrY_CallRatePerSample.txt, $!";

	print "Total of " . $ychrom . " SNP(s) available on Y chromosome.\n";

	if($ychrom > 0) {

	print "\nMendelian Error distribution on the Y chromosome based only on call rates per genotyped individual:\n\n";

		print "Pedigree\tIndividual\tSex\tProportion Y chr SNPs with called alleles\n";
		print OUT "FID\tIID\tSex\tChrY_Call_Rate\n";

		print "===========================\n";

		for($i=0;$i<$pedno;$i++){

			print "$ped[$i][0]\t$ped[$i][1]\t";
			print OUT "$ped[$i][0]\t$ped[$i][1]\t";

			if($ped[$i][4]==1){	#male

				print "male\t";
				print OUT "male\t";
			}

			if($ped[$i][4]==2){	#female

				print "female\t";
				print OUT "female\t";
			}

			if($ped[$i][4]==0){	#unknown

				print "unknown\t";
				print OUT "unknown\t";
			}

			if($whichsamples[$i]!=0){	# individual has been genotyped

				if(defined($ycheck{$ped[$i][0].$ped[$i][1]})) {

					$yprop = $ycheck{$ped[$i][0].$ped[$i][1]}/$ychromsnpcnt{$ped[$i][0].$ped[$i][1]};
					printf "%1.4f",$yprop;
					printf OUT ("%1.4f", $yprop);
				}
			}

			else {
				print "na";
				print OUT "na";
			}

			print "\n";
			print OUT "\n";

		}

	close(OUT);

	print "===========================\n";
	}

	if($males > 0 || $females > 0) {
		for($i=0;$i<$pedno;$i++){ # individual to be tested by looking for his/her parents and checking for mendelian errors
			if($whichsamples[$i]!=0){ # individual has been genotyped
				if(defined($zstat{$ped[$i][0].$ped[$i][1]})){
					if($zstat{$ped[$i][0].$ped[$i][1]}>5 || $zstat{$ped[$i][0].$ped[$i][1]} < -5){
						print_usage("Individual $ped[$i][0].$ped[$i][1] seems unlikely to have the right sex. Please check your data before proceeding.");
					}
				}
			}
		}
	}
	IND:for($i=0;$i<$pedno;$i++){ # individual to be tested by looking for his/her parents and checking for mendelian errors
#		$ind{$i}=0;
		$k=-1;
		$l=-1;
		if($whichsamples[$i]!=0){ # individual has been genotyped
			for($m=0;$m<=$#which;$m++){
				if($which[$m]==$whichsamples[$i]){
					#print "Individual $i has a corresponding which file entry of which[$m]=$which[$m]\n";
					#$m1=$m;
					if($ped[$i][2] ne "0"){ # individual is not a founder
						for($j=0;$j<$pedno;$j++){
						    if(($ped[$j][1] eq $ped[$i][2])  && ($ped[$j][0] eq $ped[$i][0])){ # find the father's genotyping sample -Tom fix
							#if($ped[$j][1] eq $ped[$i][2]){ # find the father's genotyping sample
								if($whichsamples[$j]!=0){
									#print "Found father's genotype: ped[$i][1]=$ped[$i][1],whichsamples$i]=$whichsamples[$i], ped[$i][2]=$ped[$i][2],whichsamples[$j]=$whichsamples[$j]\n";
									for($k=0;$k<=$#which;$k++){
										if($which[$k]==$whichsamples[$j]){
											#print "Individual $j (father) has a corresponding which entry of which[$k]=$which[$k]\n";
											last;
										}
									}
								}
							}
                            if(($ped[$j][1] eq $ped[$i][3])  && ($ped[$j][0] eq $ped[$i][0])){ # find the mother's genotyping sample	 - tom fix		
                            #if($ped[$j][1] eq $ped[$i][3]){ # find the mother's genotyping sample
								if($whichsamples[$j]!=0){
									#print "Found mother's genotype: ped[$i][1]=$ped[$i][1],whichsamples[$i]=$whichsamples[$i], ped[$i][3]=$ped[$i][3],whichsamples[$j]=$whichsamples[$j]\n";
									for($l=0;$l<=$#which;$l++){
										if($which[$l]==$whichsamples[$j]){
											#print "Individual $j (mother) has a corresponding which entry of which[$l]=$which[$l]\n";
											last;
										}
									}
								}
							}
						}
						if($k == -1&& $l == -1){ # test for missing father && missing mother - no test possible
							next IND;
						}
						#print "child is $m, father is $k, mother is $l\n";
						$check_poss{$i}=1;
						SNP:foreach $key (keys %genos_orig){ 
							# father's check individually
							if($k>=0){ # father
								$genof=$genos_orig{$key}[$k]+1; # note this is $k in linkdatagen also below
							}
							else{ # missing data
								$genof=0;
							}
							if($l>=0){ # mother
								$genom=$genos_orig{$key}[$l]+1;
							}
							else{ # missing data
								$genom=0;
							}
							$genoc=$genos_orig{$key}[$m]+1; 
							# work with crlmm style calls because the 0 for missing allows more compact test statements for ME checking
							#print "SNP $key: ";
							#print "genof=$genof, genom=$genom, genoc=$genoc\n";	
							# special case - X chromosome & male child
							if($annot_orig{$key}[1] eq "X" && $ped[$i][4]==1){ # X chromosome needs special checks, only carry out DAD check on autosomes, rejection ok for X chromosome, no trio check possible
								if($genom*$genoc>0){ # ensures we have no missing genotype data
									if(abs($genom-$genoc)==2){
										$ind{$i}+=1;
										$reject_snps{$key}+=1;
										$FSME{$ped[$i][0]}{$key}=1;
										next SNP;
									}
								}
							}
							else{ # do these normal checks only if either autosomes or female child
								# father's check - only autosomes
								if($genof*$genoc>0){ # ensures we have no missing genotype data
									if(abs($genof-$genoc)==2){
										$ind{$i}+=1;
										$reject_snps{$key}+=1;
										$FSME{$ped[$i][0]}{$key}=1;
										next SNP;
									}
								}
								# mother's check individually- all chromosomes
								if($genom*$genoc>0){ # ensures we have no missing genotype data
									if(abs($genom-$genoc)==2){
										$ind{$i}+=1;
										$reject_snps{$key}+=1;
										$FSME{$ped[$i][0]}{$key}=1;
										next SNP;
									}
								}
								# question: can't we still have a problem with father and son?
								# i.e. AA & BB say?
								# and a joint check
								if($genom*$genoc*$genof>0){ # ensures we have no missing genotype data
									# father AA, offspring AB, mother cannot be AA
									# father BB, offspring AB, mother cannot be BB
									if(abs($genof-$genoc)==1 && ($genof-$genom)==0){
										if($genof!=2){
											$ind{$i}+=1;
											$reject_snps{$key}+=1;
											$FSME{$ped[$i][0]}{$key}=1;
											next SNP;
										}
									}					
									# mother AA, offspring AB, father cannot be AA
									if(abs($genom-$genoc)==1 && ($genof-$genom)==0){
										if($genom!=2){
											$ind{$i}+=1;
											$reject_snps{$key}+=1;
											$FSME{$ped[$i][0]}{$key}=1;
											next SNP;
										}
									}											
								}
							}
						}
					}
				}
			}
		}
	}
	for($i=0;$i<$numchrom;$i++){
		if($i==22){
			$chr="X";
		}
		else{
			$chr=$i+1;
		}
		$chrom_tot[$i]=0;
		$chrom_mend[$i]=0;
		for($j=0;$j<$pedno;$j++){
			if($whichsamples[$j]!=0){ # individual has been genotyped
				#print "chromsnpcnt{$j}{$chr}=$chromsnpcnt{$j}{$chr}\n";
				$chromsnpcnttot{$j}+=$chromsnpcnt{$j}{$chr};
			}	
			else{
				$chromsnpcnttot{$j}=0;			
			}
		}
	}
	for($i=0;$i<=$pedno;$i++){
		$snp_cnts{$i}=0;		
	}	
	foreach $key (keys %annot_orig){
		if($annot_orig{$key}[1] eq "X"){
			$chrom=22;
			$chrom_tot[$chrom]+=1;
			if($reject_snps{$key}>0){
				$chrom_mend[$chrom]+=1;
			}
		}
		else{
			if($annot_orig{$key}[1]=~/\d/){
				$chrom=$annot_orig{$key}[1]-1;	
				$chrom_tot[$chrom]+=1;
				if($reject_snps{$key}>0){
					$chrom_mend[$chrom]+=1;
				}
			}
		}
		$snp_cnts{$reject_snps{$key}}+=1;
	}
	foreach $key (keys %annot_orig){
		if($annot_orig{$key}[1]!~/\d+/){ # not 1 -22
			if($annot_orig{$key}[1]=~/X\z/){ # X chromosome special
				$chr_snp_cnt_orig[22]+=1;
				push(@{$chr_order[22]},$key);
				#print "should be x chromo: $annot_orig{$key}[1]\n";
			}
			else{ # chr y etc and other riff raff
				$chr_snp_cnt_orig[23]+=1;
				push(@{$chr_order[23]},$key);
				#print "should be non x non auto chromo: $annot_orig{$key}[1]\n";
			}
		}
		else{
			$chr_snp_cnt_orig[$annot_orig{$key}[1]-1]+=1;
			push(@{$chr_order[$annot_orig{$key}[1]-1]},$key);
			#print "should be autosome: $annot{$key}[1]\n";
		}
	}

	print OUTM "Chr\tcM\tbp\tSNP_rsID\n";
	print OUTM "===========================\n";

	for($i=0;$i<23;$i++){

		if (%regions && !defined($regions{$i+1})){    #test if regions are being used and if so whether region included
			if ($i <= 21 || ($i == 22 && !defined($regions{"X"}))) {        ###note uses i not c
				next;
			}
		}

		%temp_chr_hash=();
 		@snporder_sorted=();

		if (defined ($chr_snp_cnt_orig[$i])) {

	 		for($j=0;$j<$chr_snp_cnt_orig[$i];$j++){
				$whichsnp=$chr_order[$i][$j];
				$whichcM=$annot_orig{$whichsnp}[7];
				$temp_chr_hash{$whichsnp}=$whichcM;
	 		# now sort into the correct order within each chromo
 			}
		}

 		sub numerically {$temp_chr_hash{$a} <=> $temp_chr_hash{$b}};
 		@snporder_sorted=sort numerically keys %temp_chr_hash;

		if (defined ($chr_snp_cnt_orig[$i])) {

			for($j=0;$j<$chr_snp_cnt_orig[$i];$j++){
				$key=$snporder_sorted[$j];

				if($reject_snps{$key}>0){

					printf OUTM ("%s\t%1.2f\t%d\t%s\n", $annot_orig{$key}[1], $annot_orig{$key}[3], $annot_orig{$key}[7], $annot_orig{$key}[0]);

#				print OUTM "$annot_orig{$key}[0] $annot_orig{$key}[1] ";
#				printf OUTM "%1.2f",$annot_orig{$key}[3]; 
#				print OUTM " ";
#				print OUTM $annot_orig{$key}[7];
#				print OUTM "\n";

				}
			}
		}
	}
	print OUTM "===========================\n";

	$chromtottot=0;
	for($i=0;$i<=22;$i++){
		$chromtottot+=$chrom_tot[$chrom];
	}
	print "\nMendelian Error distribution for genotyped individuals in the pedigree:\n\n";
	print "Individuals who have one or both parents genotyped will have Mendelian errors detected (Autosome check).\n";
	print "Individuals who are male are checked for X chromosome Mendelian errors (X chrom check).\n\n";

	open(OUT, ">$tablesDir/MendelErrorsPerSample.txt") or print_usage("Cannot write file $tablesDir/chrY_CallRatePerSample.txt, $!");

	print "Pedigree\tIndividual\tMendErrCount\tMendErrRate\tChecksConducted\n";
	print OUT "FID\tIID\tME_Freq\tME_Rate\tChecks_Conducted\n";

	print "======================\n";

	IND1:for($i=0;$i<$pedno;$i++){

		print "$ped[$i][0]\t$ped[$i][1]\t";
		print OUT "$ped[$i][0]\t$ped[$i][1]\t";

# Tests:
# 1. Is sample genotyped? If FALSE goto 4.
# 2. Is sample male? If FALSE goto 3.
#	a. Is (male) autosome NOT checked? If FALSE goto 3.
#   b. Has X chromosome been checked?
#		TRUE >> Print chrX results only.
#		FALSE >> Print "no tests possible".
# 3. Has autosome been checked?
# 		TRUE >> Print results
#		FALSE >> Print "no tests possible; likely female without parental genotype data"
# 4. >> Print "sample not genotyped".

		if($whichsamples[$i]!=0){ #test if sample genotyped

			$sum_inds+=$ind{$i};

			if($ped[$i][4]==1){   #test if male

				if(!defined($check_poss{$i})){   #test if autosome NOT checked, then can only display chrX info. if yes, then skip this.

					print "$ind{$i}\t";	#I believe this will be number of chrX mendels errors
					print OUT "$ind{$i}\t";
					$chr="X";

					if((%regions && !defined($regions{"X"})) || $chromsnpcnt{$i}{$chr} == 0) {

						print "na\tno SNPs for chromosome X available; no Mendelian error check possible\n";
						print OUT "na\tno SNPs for chromosome X available; no Mendelian error check possible\n";
					}
					
					else {					

						printf "%1.1e\tchromosome X only\n",$ind{$i}/$chromsnpcnt{$i}{$chr};
						printf OUT "%1.1e\tchromosome X only\n",$ind{$i}/$chromsnpcnt{$i}{$chr};
					}

					next IND1;
				}
			}

			if(defined($check_poss{$i})==1){	#test if autosome checked.

				print "$ind{$i}\t";		#I believe this will be number of Mendel Errors (autosome + chrX).
				print OUT "$ind{$i}\t";

				if ($chromsnpcnttot{$i}>0) {

					printf "%1.1e",$ind{$i}/$chromsnpcnttot{$i};	#I believe this will be Mendel error rate (autosome + chrX).
					printf OUT "%1.1e",$ind{$i}/$chromsnpcnttot{$i};
				}

				else {

					print "na";
					print OUT "na";
				}

				#printf "%1.2e",$ind{$i}/$chromtottot;

				print "\tautosome";
				print OUT "\tautosome";

				if($ped[$i][4]==1){

					print " and chromosome X";
					print OUT " and chromosome X";
				}

				print "\n";
				print OUT "\n";
				next IND1;
			}

			print "na\tna\tno Mendelian error check; individual likely a female without parents genotyped\n";
			print OUT "na\tna\tno Mendelian error check; individual likely a female without parents genotyped\n";
		}

		else {

			print "na\tna\tno Mendelian error check; individual not genotyped\n";
			print OUT "na\tna\tno Mendelian error check; individual not genotyped\n";
		}
	}

	close (OUT);

	print "======================\n\n";
	print "Mendelian Error distribution for genotyped SNPs:\n\n";

	print "Read as: There are <Counts> SNPs with <No of errors> (Mendelian errors).\n";
	open(OUT, ">$tablesDir/MendelErrorsPerSNP.txt") or print_usage("Cannot write file $tablesDir/MendelErrorsPerSNP.txt, $!");

	print "No of errors\tCounts\n";
	print OUT "ME_Freq\t#_SNPs\n";

	print "=======================\n";

	for($i=0;$i<=$pedno;$i++){

		print "$i\t$snp_cnts{$i}\n";
		print OUT "$i\t$snp_cnts{$i}\n";

		if($i>0){
			$mytot+=$i*$snp_cnts{$i};
		}

	}

	close(OUT);

	print "=======================\n";
	print "Total number of errors\t$mytot\n";
	print "Total number of Mendelian errors = $sum_inds\n\n";

	print "Mendelian Error distribution over all chromosomes (multiple hits per SNP not counted):\n\n";

	open(OUT, ">$tablesDir/MendelErrorsPerChromosome.txt") or print_usage("Cannot write file $tablesDir/MendelErrorsPerChromosome.txt, $!");

	print "Chr\tNo of MendErr\tMend Err Rate\n";
	print OUT "Chr\tME_Freq\tME_Rate\n";

	print "======================\n";
	for($i=0;$i<=22;$i++){

		if (%regions && !defined($regions{$i+1})){	#test if regions are being used and if so whether region included

			if ($i <= 21 || ($i == 22 && !defined($regions{"X"}))) {	###note uses i not c

				next;
			}
		}

		print $i+1,"\t$chrom_mend[$i]\t";
		print OUT $i+1 . "\t$chrom_mend[$i]\t";

		if(defined($chrom_tot[$i]) && $chrom_tot[$i]>0){

			printf "%1.4f",$chrom_mend[$i]/$chrom_tot[$i];
			printf OUT ("%1.4f",$chrom_mend[$i]/$chrom_tot[$i]);
		}
		else{
			print "0";
			print OUT "0";
		}
		print "\n";
		print OUT "\n";
	}

	close(OUT);

	print "======================\n";	

	print "Number of SNP markers prior to Mendelian error check = ",scalar keys %annot_orig,"\n";
	print "Total number of Mendelian errors detected = $mytot\n\n";


#	if (!defined($keepME)) {
	if ($MEO eq "removeSNPs") {

		print "Removing SNP markers with Mendelian errors...\n";

		foreach $key (keys %reject_snps){
			if($reject_snps{$key}>0){
				delete($genos_orig{$key});
				delete($annot_orig{$key});
				$mendelian_snp++;
			}
		}

		print "Total number of SNP markers with Mendelian errors removed = $mendelian_snp\n\n";
		print "Number of SNP markers, after removal = ",scalar keys %annot_orig,"\n";
	}
	
	elsif ($MEO eq "removeGPF") {

		print "Wiping genotypes on a per family basis...\n\n";

		open(WFSME, ">$tablesDir/familySpecificMendelianErrorGenotypesWiped.txt") or print_usage "Can't open $tablesDir/familySpecificMendelianErrorGenotypesWiped.txt, $!" ;

		print WFSME	"Family\tChr\tPhyPos\tSNP_ID\n";

		foreach $family (keys %FSME) {
		
#			print "Family: " . $family . "\n";

			foreach $key (keys %{$FSME{$family}}) {

#			print "SNP (key): " . $key . "\n";

			print "Wiping genotypes of family " . $family . " for SNP " . $key . "\n";

			print WFSME $family . "\t" . $annot_orig{$key}[1] . "\t" . $annot_orig{$key}[7] . "\t" . $key . "\n";

				foreach $i (@{$fam_i{$family}}) {

#					print "Wiping genotypes of family " . $family . " for SNP " . $key . "... individual i" . $i . "\n";
#					print "Before $genos_orig{$key}[$i] and ...";
					$genos_orig{$key}[$i] = -1 ;
#					print "... after $genos_orig{$key}[$i] \n";
				}
			}
		}
		close (WFSME);

		print "\n...done wiping genotypes per family.\n";
#		print "Total number of SNP markers with Mendelian errors removed = $mendelian_snp\n\n";
#		print "Number of SNP markers, after removal = ",scalar keys %annot_orig,"\n";
	}

	elsif ($MEO eq "keep") {

		foreach $key (keys %reject_snps){
			if($reject_snps{$key}>0){
				$mendelian_snp++;
			}
		}

		print "Total number of SNP markers with Mendelian errors = $mendelian_snp\n\n";
		print "You have chosen to keep these SNP markers (and their genotypes) with Mendelian errors.\n\n";
		print "Alternatively, you may wish to re-run LINKDATAGEN with the -MendelErrors [ removeSNPs | removeGPF ] option instead.\n";
	}

	else {
		print_usage ("Aborted due to unacceptable -MendelError value.  Allowable values are only :  keep | removeSNPs | removeGPF\n");
	}

	print "\nThere are currently:\n\n";
	print "\t", scalar keys %annot_orig, " markers from the annotation file.\n";
	print "\t", scalar keys %genos_orig, " markers from the genotype file(s).\n";
	
	close(OUTM);
}

# calculate allele frequencies and put out in a file
sub allele_freqs{
my $i;
my $key;
my $m;
my $genom;
my %allfreq=();
my $fq;
my $snpcnt=0;
my $tempsnp=0;
my $j;
my $k;
my $cnt1=0;
my $c;
my $l;
my $this;
	foreach $key(keys %annot_orig){
		$allfreq{$key}{0}=0;
		$allfreq{$key}{1}=0;
		$allfreq{$key}{2}=0;
		$allfreq{$key}{3}=0;	#to count the heteros
	}
	open(OUT, ">$tablesDir/alleleFreqs.txt") or print_usage("Can't open $tablesDir/alleleFreqs.txt, $!\n");
	IND:for($i=0;$i<$pedno;$i++){ # individual to be tested by looking for his/her parents and checking for mendelian errors
		if($whichsamples[$i]!=0){
			#print "Ind $i has genotyping data\n";
			$m=$which[$cnt1];
			for($l=0;$l<=$#which;$l++){
				if($which[$cnt1]==$which[$l]){
					$this=$l;
					last;
				}
			}
			$cnt1+=1;
			if($ped[$i][2] eq "0"){ # founder
				#print "Ind $i is a true founder\n";
				$snpcnt=0;
				foreach $key(keys %annot_orig){
					$snpcnt+=1;
					$genom=$genos_orig{$key}[$this];

					if (($ped[$i][4] == 1) && ($annot_orig{$key}[1] =~ /X\z/)) {	#special case for males and chrX

						if($genom==0){
							$allfreq{$key}{1}+=1;	#only add 1 allele count as hemizygous
						}
						if($genom==1){
#							$allfreq{$key}{1}+=1;	treat as an error and ignore
#							$allfreq{$key}{2}+=1;	treat as an error and ignore
						}
						if($genom==2){
							$allfreq{$key}{2}+=1;	#only add 1 allele count as hemizygous
						}
						if($genom== -1){
							$allfreq{$key}{0}+=1; # missing data
						}
						
					}
					
					else {

						if($genom==0){
							$allfreq{$key}{1}+=2;
						}
						if($genom==1){
							$allfreq{$key}{1}+=1;
							$allfreq{$key}{2}+=1;
							$allfreq{$key}{3}+=1;	#to count the hets
						}
						if($genom==2){
							$allfreq{$key}{2}+=2;
						}
						if($genom== -1){
							$allfreq{$key}{0}+=2; # missing data
						}
					}
				}
				next IND;
			}	
		}
	}
	print OUT "Chr\tSNP\tbp\tcM\tpopAlleleFreq\tfndrsAlleleFreq\tcntAllele1\tcntHets\tcntAllele2\tcntMiss\tpopHetero\n";
	for($c=0;$c<23;$c++){

		if (%regions && !defined($regions{$c+1})){	#test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {

				next;
			}
		}

		if (defined ($chr_snp_cnt_orig[$c])) {

			for($k=0;$k<$chr_snp_cnt_orig[$c];$k++){
				print OUT $c+1,"\t$chr_order_snp_orig[$c][$k][0]";	#chr and snp name
				print OUT "\t$chr_order_snp_orig[$c][$k][7]";		#bp
				print OUT "\t$chr_order_snp_orig[$c][$k][3]";		#cM
				printf OUT "\t%1.3f",$chr_order_snp_orig[$c][$k][4];	#populationAlleleFreq
				$key=$chr_order_snp_orig[$c][$k][0];
				if(!defined($key)){
					print_usage("chrom =$c, snp no = $k\n");
				}
				if($allfreq{$key}{1}+$allfreq{$key}{2}>0){
					$fq=$allfreq{$key}{1}/($allfreq{$key}{1}+$allfreq{$key}{2});
					printf OUT "\t%1.4lf",$fq;			#foundersAlleleFreq presumably
				}
				else{
					print OUT "\tNA";
				}
				print OUT "\t$allfreq{$key}{1}\t$allfreq{$key}{3}\t$allfreq{$key}{2}\t$allfreq{$key}{0}";	#count of allele1, hetero, allele2, followed by count of missing
				printf OUT "\t%f\n",$chr_order_snp_orig[$c][$k][5];		#heterozygosity from populationAlleleFeq
			}
		}
	}	
	close(OUT);
}
#	print "Chr\tSNP\tChr(annot)\tAllele\tPos(cM)\tAllele freq\tHet\trs\tbppos\n";

# splits all snps into chromosomes and puts them in the RIGHT order, which is necessary for 
# all .pre style output files.

# current storage in annot hash array is 0 = snpname, 1 = chromosome, 2 = allele A, 3 = cM pos, 4 =all freq A, 5 = het caucasian
sub split_snps_to_chromos{
my $i;
my $j;
my $k;
my @chr_order=(); # array of arrays with SNP names as entries
my $whichsnp;
my $whichcM;
my @snporder_sorted=();
my $key;
my $line_cnt=0;
my $snptotal=0;
my $cnt=0;
	@chr_snp_cnt_orig=(); # numbers of Snps on each chromosome for ALL Snps
	@chr_order_snp_orig=(); # array (chr) of array (genetic map) of array (all annot data) of snps ordered by 1) chromosome and 2) genetic map pos
	@chr_snp_cnt=(); # numbers of Snps on each chromosome for subset of snps
	@chr_order_snp=(); # array (chr) of array (genetic map) of array (all annot data) of snps ordered by 1) chromosome and 2) genetic map pos
	# sort into chromos first
	foreach $key (keys %annot_orig){
		#print "SNP: $key, chromo: $annot_orig{$key}[1], cM: $annot_orig{$key}[3]\n";
		if(!defined($annot_orig{$key}[1])){
			print_usage("Problem with marker $key\nSNP: $key, chromo: $annot_orig{$key}[1], cM: $annot_orig{$key}[3]\n");
		}
		if($annot_orig{$key}[1]!~/\d+/){ # not 1 -22
			if($annot_orig{$key}[1]=~/X\z/){ # X chromosome special
				$chr_snp_cnt_orig[22]+=1;
				push(@{$chr_order[22]},$key);
				#print "should be x chromo: $annot_orig{$key}[1]\n";
			}
			else{ # chr y etc and other riff raff
				$chr_snp_cnt_orig[23]+=1;
				push(@{$chr_order[23]},$key);
				#print "should be non x non auto chromo: $annot_orig{$key}[1]\n";
			}
		}
		else{
			$chr_snp_cnt_orig[$annot_orig{$key}[1]-1]+=1;
			push(@{$chr_order[$annot_orig{$key}[1]-1]},$key);
			#print "should be autosome: $annot{$key}[1]\n";
		}
		$line_cnt+=1;
	}

	open(OUT, ">$tablesDir/SNP_TotalCountPerChromosome.txt") or print_usage "Cannot open file $tablesDir/SNP_TotalCountPerChromosome.txt, $!";

	print "Chr\tTotal SNP #\n";
	print OUT "Chr\t#_SNPs\n";

	print "========================\n";
	for($i=0;$i<23;$i++){

		if (%regions && !defined($regions{$i+1})){	#test if regions are being used and if so whether region included

			if ($i <= 21 || ($i == 22 && !defined($regions{"X"}))) {	###note uses i not c

				next;
			}
		}

		if (defined ($chr_snp_cnt_orig[$i])) {

			print $i+1,"\t$chr_snp_cnt_orig[$i]\n";
			print OUT $i+1 . "\t$chr_snp_cnt_orig[$i]\n";

			$snptotal+=$chr_snp_cnt_orig[$i];
		}

		else {

			print $i+1, "\t0\n";
			print OUT $i+1, "\t0\n";
		}

	}
	close(OUT);

	print "========================\n";
	print "Total\t$snptotal (not including chrY SNPs)\n";
	# sorting the chromosome specific arrays by genetic map position   #CATHEDIT: sort by bp pos

	print "\nSorting SNPs by genetic position ...\n";

	for($i=0;$i<23;$i++){

		if (%regions && !defined($regions{$i+1})) {

			if ($i <= 21 || ( $i == 22 && !defined($regions{"X"}))){

				next;
			}
		}

		print 
		%temp_chr_hash=();
		@snporder_sorted=();

		if (defined ($chr_snp_cnt_orig[$i])) {

			for($j=0;$j<$chr_snp_cnt_orig[$i];$j++){
				$whichsnp=$chr_order[$i][$j];
				$whichcM=$annot_orig{$whichsnp}[7];  #CATHEDIT: used to be {$whichsnp}[3]
				$temp_chr_hash{$whichsnp}=$whichcM;
			# now sort into the correct order within each chromo
			}
			#sub numerically {$temp_chr_hash{$a} <=> $temp_chr_hash{$b}};
			@snporder_sorted=sort numerically keys %temp_chr_hash;
			# creating the final, sorted array (chr) of arrays (snps) of arrays (all the data we need)
			for($j=0;$j<$chr_snp_cnt_orig[$i];$j++){
				$key=$snporder_sorted[$j];
				for($k=0;$k<=7;$k++){
					$chr_order_snp_orig[$i][$j][$k]=$annot_orig{$key}[$k];
				}
			}
		}
	}
}

sub print_ordered_snps{
my $i;
my $j;
my $k;
	open(OUT10, ">$tablesDir/orderedSNPs.txt") or print_usage "Cannot open file $tablesDir/orderedSNPs.txt, $!";
	print OUT10 "Chr\tSNP\tChr(annot)\tAllele\tPos(cM)\tAllele freq\tHet\trs\tbppos\n";
	for($i=0;$i<23;$i++){

		if (%regions && !defined($regions{$i+1})){	#test if regions are being used and if so whether region included

			if ($i <= 21 || ($i == 22 && !defined($regions{"X"}))) {	###note uses i not c

				next;
			}
		}

		if (defined ($chr_snp_cnt_orig[$i])) {

			for($j=0;$j<$chr_snp_cnt_orig[$i];$j++){
				print OUT10 $i+1,"\t";
				for($k=0;$k<=7;$k++){
					print OUT10 "$chr_order_snp_orig[$i][$j][$k]\t";
				}
				print OUT10 "\n";
			}

			print OUT10 "-------------------------\n";
			print OUT10 "$chr_snp_cnt_orig[$i]\n";
			print OUT10 "-------------------------\n";
		}
	}
}

# --------subsetting of the data --------------------------------------

# function to select only a subset of markers 
# allow for two options: 1) random choice of marker per bin or 2) pick the "best" in the sense of max het.
sub select_per_bin{

%genos=();

my $i;
my $j;
my $k;
my $entry;
my $key;
my $value;
my $rand_snp;
my @temp_snp_order=();
my @final_list_of_snps=();
my @final_list_of_snps_names=();
my @current_bin_markers=();
my $snp_length=0;
my @snporder_sorted=();
my $current_marker_start=0;
my @keyorder=();
my $sum_of_snps=();
my $sum_of_snps_orig=();
my @final_list_of_hets=();
my $avge_hets=0;
my $mytest=1; # for error checking within this function
#CATHEDIT: This parameter $mytest seems redundant.  Setting it to 0 causes downstream errors.
#Quick patch is to set this to 1, this used to be $mytest=0.
my $assigned; #CATHEDIT: for preventing SNPs from being added to the list if they are less than minDist
				#away from previous SNP
my $mydist;
my $chromosome;
my $newmappos;
my $oldmappos;
my $test=0;
my $final=0;
my $totalHetero=0;
my $totalCountSNPs=0;
my $randOffSet = -1;
my $randRange;

	# subselect in chr_order_snp & reduce chr_snp_cnt too (these contain all the necessary information
	print OUT1 "Chromosome\tInterval\n";
	for($i=0;$i<23;$i++){

		if ( (%regions && !defined($regions{$i+1})) || !defined($chr_snp_cnt_orig[$i]) ){	#test if regions are being used and if so whether region included

			if ($i <= 21 || ($i == 22 && !defined($regions{"X"})) || !defined($chr_snp_cnt_orig[$i]) ) {	###note uses i not c

				next;
			}
		}

		for($j=0;$j<$chr_snp_cnt_orig[$i];$j++){
			for($k=0;$k<=$#which;$k++){
				#print "assigning step, k = ".$k."\n";
				if(!defined($genos_orig{$chr_order_snp_orig[$i][$j][0]}[$k])){  # SPECIAL CLAUSE FOR MY SiMULATIONS
					$genos_orig{$chr_order_snp_orig[$i][$j][0]}[$k]=0;
				}
			}
		}
	}
	if($binsize==0){ # special case, don't select any snps.
		for($i=0;$i<23;$i++){

			if ((%regions && !defined($regions{$i+1})) || !defined($chr_snp_cnt_orig[$i]) ){	#test if regions are being used and if so whether region included

				if ($i <= 21 || ($i == 22 && !defined($regions{"X"}))) {	###note uses i not c

					next;
				}
			}

			$chr_snp_cnt[$i]=$chr_snp_cnt_orig[$i];

			for($j=0;$j<$chr_snp_cnt_orig[$i];$j++){
				for($k=0;$k<=7;$k++){
					$chr_order_snp[$i][$j][$k]=$chr_order_snp_orig[$i][$j][$k];
				}
				for($k=0;$k<=$#which;$k++){
					$genos{$chr_order_snp[$i][$j][0]}[$k]=$genos_orig{$chr_order_snp[$i][$j][0]}[$k];
				}
			}
		}

		return();
	}

	print "\nDistribution of SNPs selected per chromosome:\n\n";

	open(OUT, '>', "$tablesDir/SNP_SelectedCountPerChromosome.txt") or print_usage "Cannot open file $tablesDir/SNP_SelectedCountPerChromosome.txt, $!";

	print "Chr\t# Selected SNPs\tAverage Heterozygosity\n";
	print OUT "Chr\t#_SNPs\tAverage_Heterozygosity\n";

	print "========================\n";

	for($i=0;$i<23;$i++){

		if (%regions && !defined($regions{$i+1})){	#test if regions are being used and if so whether region included

			if ($i <= 21 || ($i == 22 && !defined($regions{"X"}))) {	###note uses i not c and

				next;
			}
		}

		$avge_hets=0;
		if($mytest==0){
			print "Chromosome ",$i+1,"\n";
		}
		$snp_length=$binsize;
		$current_marker_start=0;
		$j=0;

		if (defined ($chr_snp_cnt_orig[$i])) {

		SECTION: while($j<$chr_snp_cnt_orig[$i]){			#each loop of this SECTION represents a new bin
			# find the number of SNPs in the bin of length x cM
			if($chr_order_snp_orig[$i][$j][3]<$snp_length){
				if(defined($randomSNP)){
					push(@current_bin_markers,$j);
					#print "current_bin_markers=@current_bin_markers\n";
				}
				else{
					$bin_het_hash{$j}=$chr_order_snp_orig[$i][$j][5];
				}
				$j+=1;
				if($j<($chr_snp_cnt_orig[$i]-1)){
					next SECTION; # if last marker continue to below, else keep stacking markers on.
				}
			}
			# end of interval sized bin_size: assess which SNP to pick, either random or highest het.
			$snp_length+=$binsize;
			# option 1: random snp:
			if(defined($randomSNP) && $mytest==1){
				if($#current_bin_markers == -1){
					print OUT1 $i+1,"No snp in interval (",$snp_length-$binsize,",",$snp_length,")\n";
					next SECTION;
				}
				$final=0;
					$assigned = 0;

					## loop through the SNPs until find SNP that is greater in cM than minimum dist from last selected SNP (oldpos).
					## then choose random SNP from set of SNPs including that SNP plus all remaining SNPs in bin

					if($#final_list_of_snps == -1){
						$oldmappos=0;
					}
					else{
						$oldmappos=$chr_order_snp_orig[$i][$final_list_of_snps[$#final_list_of_snps]][3];
					}

					for($k=0;$k<=$#current_bin_markers;$k++){	#loop through SNPs in bin

						if ( ($chr_order_snp_orig[$i][$current_bin_markers[$k]][3] - $oldmappos) >= $minDist ) {

							$randOffSet = $k;	#find SNP that is closest to last SNP selected, but greater than minDist; 
							last;
						}
					}
					
					if ($randOffSet == -1) {
					
##TEMP						print "There is no marker in this bin with distance from previous SNP less than $minDist .  Last marker dist is $chr_order_snp_orig[$i][$#current_bin_markers+$current_bin_markers[0]][3]\n";
						@current_bin_markers=();
						next SECTION;	#goto next bin
					}
					else {

						$randRange = ((scalar @current_bin_markers) - $randOffSet);
					}

					$rand_snp = rand($randRange); # chooses a random number from the range of suitable SNPs in the bin
					$rand_snp = int($rand_snp); # must make into integer
					$final = $rand_snp + $randOffSet;
					$assigned = 1;

				if($assigned == 1) {
					push(@final_list_of_snps,$final+$current_bin_markers[0]); # rescale to snp no in chromo
					push(@final_list_of_snps_names,$chr_order_snp_orig[$i][$final+$current_bin_markers[0]][0]);
					push(@final_list_of_hets,$chr_order_snp_orig[$i][$final+$current_bin_markers[0]][5]); # snp het.. 
					}
				if($mytest==0){
					print "SNP choice is random:\n";
					for($k=0;$k<=$#current_bin_markers;$k++){
						print "$current_bin_markers[$k] ";
					}
					print "\n";
				}
				
				$randOffSet = -1;
				
				@current_bin_markers=();
			}
			# option 2: choose SNP with highest het (according to CEPH (or other) HAPMAP data)
			else{ 
				if(scalar keys(%bin_het_hash) == 0){
					printf OUT1 $i+1,"No snp in interval (",$snp_length-$binsize,",",$snp_length,")\n";
					next SECTION;
				}
				if($mytest==0){
					print "SNP choice is highest het:\n";
					while (($key,$value) = each (%bin_het_hash)){
						printf "%1.2f",$key;
						print ", ";
						printf "%1.2f",$value;
						print "\n";
					}
				}
				$current_marker_start=$j;
				sub binnumerically {$bin_het_hash{$b} <=> $bin_het_hash{$a}}; # b and a reversed: sort is in descending order
				@snporder_sorted=sort binnumerically keys %bin_het_hash;
				$chromosome = $i+1;
				print OUT1 "chr".$chromosome." bin (",$snp_length-$binsize,",",$snp_length,")\n";
				#print $chr_order_snp_orig[$i][$snporder_sorted[$final]][0]."\n";
				#print "highest het snp is $snporder_sorted[0], het = $bin_het_hash{$snporder_sorted[0]}\n";
				# now making sure that the mrkers have a minimum distance apart, at least 0.2 cM
				$final=0;
				if($minDist>0){
					#print "oldmappos = ".$oldmappos."\n";
					$assigned = 0;
					for($k=0;$k<=$#snporder_sorted;$k++){
#					for($k=$#snporder_sorted;$k>=0;$k--){	# if wish to select the least het
						#print "k = ".$k."\n";
						$newmappos=$chr_order_snp_orig[$i][$snporder_sorted[$k]][3];
						#print "In minDist loop: newmapppos = ".$newmappos."\n";
						if($#final_list_of_snps == -1){
							$oldmappos=0;
						}
						else{
							$oldmappos=$chr_order_snp_orig[$i][$final_list_of_snps[$#final_list_of_snps]][3];
						}
						if(!defined($chr_order_snp_orig[$i][$snporder_sorted[$k]][3])){
							print_usage("chr_order_snp_orig[$i][$final_list_of_snps[$#final_list_of_snps]][3]=$chr_order_snp_orig[$i][$final_list_of_snps[$#final_list_of_snps]][3]\n");
						}
						$mydist=$newmappos-$oldmappos;
						#print "old ".$oldmappos."   new ".$newmappos." dist ".$mydist."\n";
						if($mydist >= $minDist){
							#print "finalloop = ".$k."\n";
							$final=$k;
							$assigned = 1;
							last;
						}
					}
				}
				#print "\#snpordersorted ".$#snporder_sorted."\n";
				#if($#snporder_sorted==0){
					if($mydist<$minDist){
						printf OUT1 "There is no marker in this bin with distance from previous SNP greater than or equal to $minDist cM.\n";
						#printf "%1.2f",$mydist;
						#print "cM.\n";
					#}
				}
				if($assigned == 1) {
					push(@final_list_of_snps,$snporder_sorted[$final]); # choose HIghest het snp 
					#print "final = ".$final."\n";
					#print "push ".$chr_order_snp_orig[$i][$snporder_sorted[$final]][0]."\n";
					push(@final_list_of_snps_names,$chr_order_snp_orig[$i][$snporder_sorted[$final]][0]); # snp name.. 

					#print "chr $i final $final snp_name " . $chr_order_snp_orig[$i][$snporder_sorted[$final]][0] . " ";
					#print "heterozygosity: " . $chr_order_snp_orig[$i][$snporder_sorted[$final]][5] . "\n";

					push(@final_list_of_hets,$chr_order_snp_orig[$i][$snporder_sorted[$final]][5]); # snp het.. 
				}
			}
			%bin_het_hash=();
		}				# END OF SECTION curly bracket - next bin or carry on...

		}

		if($mytest==0){
			print "Chromosome ",$i+1,":\n";
			foreach $entry (@final_list_of_snps){
				print "$entry\t";
				for($k=0;$k<=7;$k++){
					if(defined($chr_order_snp_orig[$i][$entry][$k])){
						print "$chr_order_snp_orig[$i][$entry][$k]\t";
					}
					else{
						print "\nProblem with [$i][$entry][$k], total count on chromosome is $chr_snp_cnt_orig[$i]\n";
					}
				}
				print "\n";
			}
		}
		# now need to only keep the markers listed in final_list_of_snps and ditch the rest
		# first a check however.. print list of snps to use:

		print $i+1 . "\t" . ($#final_list_of_snps + 1) . "\t";

		print OUT $i+1 . "\t" . ($#final_list_of_snps + 1) . "\t";

		print OUT1 "# of Snps chosen from chromosome ",$i+1," is ",$#final_list_of_snps+1,"\n";

		# writing out the new data:
		$chr_snp_cnt[$i]=$#final_list_of_snps+1;
		$sum_of_snps+=$chr_snp_cnt[$i];

		if (defined($chr_snp_cnt_orig[$i])) {

			$sum_of_snps_orig+=$chr_snp_cnt_orig[$i];
		}

		for($j=0;$j<=$#final_list_of_snps;$j++){
			$avge_hets+=$final_list_of_hets[$j];
			$totalHetero+=$final_list_of_hets[$j];
			for($k=0;$k<=7;$k++){
				$chr_order_snp[$i][$j][$k]=$chr_order_snp_orig[$i][$final_list_of_snps[$j]][$k];
			}
			for($k=0;$k<=$#which;$k++){
				$genos{$final_list_of_snps_names[$j]}[$k]=$genos_orig{$final_list_of_snps_names[$j]}[$k];
			}
		}

		if($#final_list_of_snps != -1){

			$totalCountSNPs +=scalar(@final_list_of_snps);
			$avge_hets/=($#final_list_of_snps+1);
			printf "%1.4f", $avge_hets;
			print "\n";
			printf OUT ("%1.4f\n", $avge_hets);

			print OUT1 "Average heterozygosity of chosen snps = ";
			printf OUT1 "%1.4f", $avge_hets;
			print OUT1 "\n\n";
		}
		else {
			print "na\n";
			print OUT "na\n";
			print OUT1 "Average heterozygosity of chosen snps = na\n\n";
		}

		@final_list_of_snps=();
		@final_list_of_snps_names=();
		@final_list_of_hets=();
	}

	close(OUT);

	print "======================\n";

	if (defined($sum_of_snps_orig)) {

		print "Number of SNPs selected: $sum_of_snps, from a total number of original snps: $sum_of_snps_orig\n";
		printf ("\nAverage heterozygosity of these selected SNPs is: %1.4f\n", ($totalHetero/$totalCountSNPs));
	}

	else {

		print "No SNPs chosen as no original SNPs exist to choose from.\n\n";
	}

		foreach $key (keys %genos_orig) {
	
			if( !defined($genos{$key}) ){
		
				delete($genos_orig{$key});
				delete($annot_orig{$key}); 
			}
		}

	print "\nThere are currently:\n\n";
	print "\t", scalar keys %annot_orig, " markers from the annotation file.\n";
	print "\t", scalar keys %genos_orig, " markers from the genotype file(s).\n";
}

sub popHetTest {

	#compare the homozygous/heterozygous genotype counts of SNPs against the 'selected' population estimates

	my $i;
	my $j;
	my $chr;
	my @sex = ();
	my $popHomoAA = 0;
	my $popHomoBB = 0;
	my $popHeteroAB = 0;
	my $geno = 0;
	my @chiSq = (); #not set to zero yet
	my @chiSqAutosome = (); #not set to zero yet
	my @chiSqGenome = (); #not set to zero yet
	my @chrCountSNPs = ();
	my $genomeCountSNPs = 0;
	my $autosomeCountSNPs = 0;

	##for population based counts and sums:

	my @popChrSumHomoAA = ();
	my @popChrSumHomoBB = ();
	my @popChrSumHeteroAB = ();
	my @popChrMeanHet = ();

	my $popGenomeSumHomoAA = 0;		# is also the expected pop HomoAA
	my $popGenomeSumHomoBB = 0;		# is also the expected pop HomoBB
	my $popGenomeSumHeteroAB = 0;	# is also the expected pop HeteroAB
	my $popGenomeMeanHet = 0;

	##for individual specific counts and sums:

	my @indChrGenoCount = ();
	my @indChrMissCount = ();
	my @indChrCountHomoAA = ();
	my @indChrCountHomoBB = ();
	my @indChrCountHeteroAB = ();
	
	my @indGenomeGenoCount = ();
	my @indGenomeMissCount = ();
	my @indGenomeCountHomoAA = (); 
	my @indGenomeCountHomoBB = ();
	my @indGenomeCountHeteroAB = ();
	
	my @indPopChrSumHomoAA = ();
	my @indPopChrSumHomoBB = ();
	my @indPopChrSumHeteroAB = ();

	my @indPopGenomeSumHomoAA = ();
	my @indPopGenomeSumHomoBB = ();
	my @indPopGenomeSumHeteroAB = ();

	##for autosome only, population:

	my $popAutosomeExpHomoAA = 0;
	my $popAutosomeExpHomoBB = 0;
	my $popAutosomeExpHeteroAB = 0;
	my $popAutosomeMeanHet = 0;

	##for autosome only, individual specific:

	my $indAutosomeGenoCount = 0;
	my $indAutosomeMissCount = 0;
	my $indAutosomeCountHomoAA = 0;
	my $indAutosomeCountHomoBB = 0;
	my $indAutosomeCountHeteroAB = 0;

	my @indAutosomeExpHomoAA = (); #not set to zero yet
	my @indAutosomeExpHomoBB = (); #not set to zero yet
	my @indAutosomeExpHeteroAB = (); #not set to zero yet

	if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {
	
		open(OUT, ">$tablesDir/popHetTestDetails.txt") or print_usage "Cannot open file $tablesDir/popHetTestDetails.txt, $!";
	}

	print "\nCommencing goodness-of-fit test for observed homozygous (AA & BB) and heterozygous (AB) genotype counts given the selected population allele frequencies ...\n";

	#### printing header information for table - length is dependent on number of genotyped individuals

	if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

		print "\n\t\t\t\t\t";
		print OUT "\n\t\t\t\t\t";
	}

	for ($i=0; $i<$pedno; $i++){

		for($j=0;$j<=$#which;$j++){ # genotyped individual

			if($which[$j]==$whichsamples[$i]){

				if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

					if ( length($ped[$i][0]) > 7 ) {

						print "\t|----\t\tFamily:\t" . substr($ped[$i][0], 0, 7) . "...\t  ----|";
					}

					else {
				
						print "\t|----\t\tFamily:\t$ped[$i][0]\t\t  ----|";
					
					}

					print OUT "\t|----\t\tFamily:\t$ped[$i][0]\t\t  ----|";
				}

				##setting individual genome counts to zero
				$indGenomeGenoCount[$i] = 0;
				$indGenomeMissCount[$i] = 0;
				$indGenomeCountHomoAA[$i] = 0;
				$indGenomeCountHomoBB[$i] = 0;
				$indGenomeCountHeteroAB[$i] = 0;
				$indPopGenomeSumHomoAA[$i] = 0;
				$indPopGenomeSumHomoBB[$i] = 0;
				$indPopGenomeSumHeteroAB[$i] = 0;
			}
		}
	}

	if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

		print "\n\t\t\t\t\t";
		print OUT "\n\t\t\t\t\t";
	}

	for ($i=0; $i<$pedno; $i++){

		for($j=0;$j<=$#which;$j++){ # genotyped individual

			if($which[$j]==$whichsamples[$i]){
			
				if ($ped[$i][4] == 1) {
					$sex[$i] = "Male";
				}
				elsif ($ped[$i][4] == 2) {
					$sex[$i] = "Female";
				}
				else {
					$sex[$i] = "Unkown";
				}

				if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

					if ( length($ped[$i][1]) > 7 ) {		#test for header length. if too long then needs one less tab for onscreen.

						print "\t|----\tMember:\t" . substr($ped[$i][1], 0, 7) . "... sex:\t$sex[$i]\t  ----|";
					}

					else {

						print "\t|----\tMember:\t$ped[$i][1]\t   sex:\t$sex[$i]\t  ----|";
					}
				
					print OUT "\t|----\tMember:\t$ped[$i][1]\tsex:\t$sex[$i]\t  ----|";
				}
			}
		}
	}

	if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

		print "\nChr(s)\t #SNPs\t$pop" . "_Het\tExHomAA\tExHetAB\tExHomBB";
		print OUT "\nChromosome\t#SNPs\t$pop" . "_Heterozygosity\tExHomAA\tExHetAB\tExHomBB";
	}

	for ($i=0; $i<$pedno; $i++){

		for($j=0;$j<=$#which;$j++){ # genotyped individual 

			if($which[$j]==$whichsamples[$i]){

				if ($popHetTest eq "verbose" || $popHetTest eq "perChrVerbose") {

					print "\tMissing\tGeno\tHom_AA\tHet_AB\tHom_BB\t chi^2";
					print OUT "\tMissing\tGeno\tHom_AA\tHet_AB\tHom_BB\t chi^2";
				}

				elsif ($popHetTest eq "perChr") {

					print "\tObsGeno\tObsMiss\tObHomAA\tObHetAB\tObHomBB\t chi^2";
					print OUT "\tObsGeno\tObsMiss\tObHomAA\tObHetAB\tObHomBB\t chi^2";
				}
			}
		}
	}

	if ($popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

		print "\n-------\t-------\t-------\t-------\t-------\t-------";
		print OUT "\n-------\t-------\t-------\t-------\t-------\t-------";
	
		for ($i=0; $i<$pedno; $i++){

			for($j=0;$j<=$#which;$j++){ # genotyped individual

				if($which[$j]==$whichsamples[$i]){

					print "\t-------\t-------\t-------\t-------\t-------\t-------";
					print OUT "\t-------\t-------\t-------\t-------\t-------\t-------";
				}
			}
		}
	}

	#### now working out the expected and observed homozygosity and heterozygosity rates

	for ($chr = 1; $chr <= $numchrom; $chr++) {		#loop to set things to zero

		$chrCountSNPs[$chr] = 0;
		$popChrSumHomoAA[$chr] = 0;
		$popChrSumHomoBB[$chr] = 0;
		$popChrSumHeteroAB[$chr] = 0;

		for ($i=0; $i<$pedno; $i++){

			for ($j=0; $j<=$#which; $j++){	

				if($which[$j] == $whichsamples[$i]){

					$indPopChrSumHomoAA[$chr][$i] = 0;
					$indPopChrSumHomoBB[$chr][$i] = 0;
					$indPopChrSumHeteroAB[$chr][$i] = 0;
					$indChrGenoCount[$chr][$i] = 0;
					$indChrMissCount[$chr][$i] = 0;
					$indChrCountHomoAA[$chr][$i] = 0;
					$indChrCountHomoBB[$chr][$i] = 0;
					$indChrCountHeteroAB[$chr][$i] = 0;
				}
			}
		}
	}

	foreach $key (keys %annot_orig) {	#second loop to add up the heterozygosity, and heterozygous and homozygous genotypes

		$chr = $annot_orig{$key}[1];	

		if ($chr eq "X") {
		
			$chr = 23;
		}
		elsif ($chr eq "XY" || $chr eq "Y" || $chr eq "M") {

			next;
		}
		elsif ($chr =~ /\D/) {

			die("Problem with annotation file, this value should be a chromosome number: $chr");
		}

		++$chrCountSNPs[$chr];
		++$genomeCountSNPs;
		
		$popHomoAA = $annot_orig{$key}[4] ** 2;
		$popHomoBB = (1 - $annot_orig{$key}[4]) ** 2;
		$popHeteroAB = 1 - ($popHomoAA + $popHomoBB);

		$popChrSumHomoAA[$chr] += $popHomoAA;
		$popChrSumHomoBB[$chr] += $popHomoBB;
		$popChrSumHeteroAB[$chr] += $popHeteroAB;

		$popGenomeSumHomoAA += $popHomoAA;
		$popGenomeSumHomoBB += $popHomoBB;
		$popGenomeSumHeteroAB += $popHeteroAB;
		
		for ($i=0; $i<$pedno; $i++){

			for ($j=0; $j<=$#which; $j++){	

				if($which[$j] == $whichsamples[$i]){
				
					$geno = $genos_orig{$key}[$j];
						
					if ($geno == -1) {

						++$indChrMissCount[$chr][$i];
						++$indGenomeMissCount[$i];
					}
					
					else {

						++$indChrGenoCount[$chr][$i];
						++$indGenomeGenoCount[$i];

						if ($ped[$i][4] == 1 && $chr == 23) {	### if male and chr==23 then do special thing, otherwise do normal

							if ($geno == 1) {

								++$indChrMissCount[23][$i];	#chr = 23
								++$indGenomeMissCount[$i];

								$indPopChrSumHeteroAB[$chr][$i] = "NA:male";	#I know this is silly, as it occurs for every SNP on chrX, but is there for emphasis;
							}

							else {

								$indPopChrSumHomoAA[23][$i] += $annot_orig{$key}[4];	#chr = 23
								$indPopGenomeSumHomoAA[$i] += $annot_orig{$key}[4];

								$indPopChrSumHomoBB[23][$i] += (1 - $annot_orig{$key}[4]);	#chr = 23
								$indPopGenomeSumHomoBB[$i] += (1 - $annot_orig{$key}[4]);
							}
						}

						else {

							$indPopChrSumHomoAA[$chr][$i] += $popHomoAA;
							$indPopGenomeSumHomoAA[$i] += $popHomoAA;

							$indPopChrSumHeteroAB[$chr][$i] += $popHeteroAB;
							$indPopGenomeSumHeteroAB[$i] += $popHeteroAB;
						
							$indPopChrSumHomoBB[$chr][$i] += $popHomoBB;
							$indPopGenomeSumHomoBB[$i] += $popHomoBB;
						}

						if ($geno == 0) {

							++$indChrCountHomoAA[$chr][$i];
							++$indGenomeCountHomoAA[$i];
						}

						if ($geno == 1) {

							++$indChrCountHeteroAB[$chr][$i];
							
							if (!($ped[$i][4] == 1 && $chr == 23)) {  # if not male and chr23

								++$indGenomeCountHeteroAB[$i];
							}
						}

						if ($geno == 2) {

							++$indChrCountHomoBB[$chr][$i];
							++$indGenomeCountHomoBB[$i];
						}
					}
				}
			}
		}
	}
	
	for($chr = 1; $chr <= $numchrom; $chr++){

		if ( (%regions && !defined($regions{$chr})) ) {   ##|| !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($chr <= 22 || ($chr == 23 && !defined($regions{"X"}))) {
				next;
			}
		}

		if ($chrCountSNPs[$chr] > 0) {

			$popChrMeanHet[$chr] = $popChrSumHeteroAB[$chr] / $chrCountSNPs[$chr];
		}
		
		else {
		
			$popChrMeanHet[$chr] = -9;
		}

		if ($popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

			printf ("\n%4d\t%7d\t%7.4f\t%7.1f\t%7.1f\t%7.1f", $chr, $chrCountSNPs[$chr], $popChrMeanHet[$chr], $popChrSumHomoAA[$chr], $popChrSumHeteroAB[$chr], $popChrSumHomoBB[$chr]);
			printf OUT ("\n%4d\t%7d\t%7.4f\t%7.1f\t%7.1f\t%7.1f", $chr, $chrCountSNPs[$chr], $popChrMeanHet[$chr], $popChrSumHomoAA[$chr], $popChrSumHeteroAB[$chr], $popChrSumHomoBB[$chr]);
		}

		for ($i=0; $i<$pedno; $i++) {

			for($j=0; $j<=$#which; $j++){ # genotyped individual

				if($which[$j]==$whichsamples[$i]){

					if ($ped[$i][4] == 1 && $chr == 23) {  #if male and chr23 then ignore the heteroAB counts

						if (($indPopChrSumHomoAA[$chr][$i] * $indPopChrSumHomoBB[$chr][$i]) > 0) {

							$chiSq[$i] = (($indChrCountHomoAA[$chr][$i] - $indPopChrSumHomoAA[$chr][$i]) ** 2 / $indPopChrSumHomoAA[$chr][$i]) + (($indChrCountHomoBB[$chr][$i] - $indPopChrSumHomoBB[$chr][$i]) ** 2 / $indPopChrSumHomoBB[$chr][$i]);
						}
						
						else {

							$chiSq[$i] = -9;
						}
					}
					
					else {

						if (($indPopChrSumHomoAA[$chr][$i] * $indPopChrSumHeteroAB[$chr][$i] * $indPopChrSumHomoBB[$chr][$i]) > 0) {

							$chiSq[$i] = (($indChrCountHomoAA[$chr][$i] - $indPopChrSumHomoAA[$chr][$i]) ** 2 / $indPopChrSumHomoAA[$chr][$i]) + (($indChrCountHeteroAB[$chr][$i] - $indPopChrSumHeteroAB[$chr][$i]) ** 2 / $indPopChrSumHeteroAB[$chr][$i]) + (($indChrCountHomoBB[$chr][$i] - $indPopChrSumHomoBB[$chr][$i]) ** 2 / $indPopChrSumHomoBB[$chr][$i]);
						}
						
						else {

							$chiSq[$i] = -9;
						}
					}

					if ($popHetTest eq "perChrVerbose") {

						printf ("\t%7d\tObserve\t%7d\t%7d\t%7d\t", $indChrMissCount[$chr][$i], $indChrCountHomoAA[$chr][$i], $indChrCountHeteroAB[$chr][$i], $indChrCountHomoBB[$chr][$i]);
						printf OUT ("\t%7d\tObserve\t%7d\t%7d\t%7d\t", $indChrMissCount[$chr][$i], $indChrCountHomoAA[$chr][$i], $indChrCountHeteroAB[$chr][$i], $indChrCountHomoBB[$chr][$i]);
					}

					elsif ($popHetTest eq "perChr") {

						printf ("\t%7d\t%7d\t%7d\t%7d\t%7d\t%7.2f", $indChrGenoCount[$chr][$i], $indChrMissCount[$chr][$i], $indChrCountHomoAA[$chr][$i], $indChrCountHeteroAB[$chr][$i], $indChrCountHomoBB[$chr][$i], $chiSq[$i]);
						printf OUT ("\t%7d\t%7d\t%7d\t%7d\t%7d\t%7.2f", $indChrGenoCount[$chr][$i], $indChrMissCount[$chr][$i], $indChrCountHomoAA[$chr][$i], $indChrCountHeteroAB[$chr][$i], $indChrCountHomoBB[$chr][$i], $chiSq[$i]);
					}
				}
			}
		}

		if ($popHetTest eq "perChrVerbose") {
		
			print "\n\t\t\t\t\t";
			print OUT "\n\t\t\t\t\t";

			for ($i=0; $i<$pedno; $i++) {

				for($j=0; $j<=$#which; $j++){ # genotyped individual

					if($which[$j]==$whichsamples[$i]){

						if ($ped[$i][4] == 1 && $chr == 23) {  #if male and chr23 then print differently

							printf ("\t\tExpect\t%7.1f\t%7s\t%7.1f\t%7.2f", $indPopChrSumHomoAA[$chr][$i], $indPopChrSumHeteroAB[$chr][$i], $indPopChrSumHomoBB[$chr][$i], $chiSq[$i]);
							printf OUT ("\t\tExpect\t%7.1f\t%7s\t%7.1f\t%7.2f", $indPopChrSumHomoAA[$chr][$i], $indPopChrSumHeteroAB[$chr][$i], $indPopChrSumHomoBB[$chr][$i], $chiSq[$i]);
						}
						
						else {

							printf ("\t\tExpect\t%7.1f\t%7.1f\t%7.1f\t%7.2f", $indPopChrSumHomoAA[$chr][$i], $indPopChrSumHeteroAB[$chr][$i], $indPopChrSumHomoBB[$chr][$i], $chiSq[$i]);
							printf OUT ("\t\tExpect\t%7.1f\t%7.1f\t%7.1f\t%7.2f", $indPopChrSumHomoAA[$chr][$i], $indPopChrSumHeteroAB[$chr][$i], $indPopChrSumHomoBB[$chr][$i], $chiSq[$i]);
						}
					}
				}
			}					
		}
	}

	$autosomeCountSNPs = $genomeCountSNPs - $chrCountSNPs[23];

	$autosomeCountSNPs{$popCol} = $autosomeCountSNPs;

	if ($autosomeCountSNPs > 0) {

		$popAutosomeExpHomoAA = $popGenomeSumHomoAA - $popChrSumHomoAA[23];
		$popAutosomeExpHomoBB = $popGenomeSumHomoBB - $popChrSumHomoBB[23];
		$popAutosomeExpHeteroAB = $popGenomeSumHeteroAB - $popChrSumHeteroAB[23];
		$popAutosomeMeanHet = $popAutosomeExpHeteroAB / $autosomeCountSNPs;

		$popAutosomeMeanHet{$popCol} = $popAutosomeMeanHet;

		if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

			print "\n-------\t-------\t-------\t-------\t-------\t-------";
			print OUT "\n-------\t-------\t-------\t-------\t-------\t-------";

			for ($i=0; $i<$pedno; $i++){

				for($j=0;$j<=$#which;$j++){ # genotyped individual

					if($which[$j]==$whichsamples[$i]){

						print "\t-------\t-------\t-------\t-------\t-------\t-------";
						print OUT "\t-------\t-------\t-------\t-------\t-------\t-------";
					}
				}
			}

			printf ("\nautosme\t%7d\t%7.4f\t%7.1f\t%7.1f\t%7.1f", $autosomeCountSNPs, $popAutosomeMeanHet, $popAutosomeExpHomoAA, $popAutosomeExpHeteroAB, $popAutosomeExpHomoBB);
			printf OUT ("\nautosome\t%7d\t%7.4f\t%7.1f\t%7.1f\t%7.1f", $autosomeCountSNPs, $popAutosomeMeanHet, $popAutosomeExpHomoAA, $popAutosomeExpHeteroAB, $popAutosomeExpHomoBB);
		}

		for ($i=0; $i<$pedno; $i++){

			for ($j=0; $j<=$#which; $j++){	

				if($which[$j] == $whichsamples[$i]){

					$indAutosomeGenoCount = $indGenomeGenoCount[$i] - $indChrGenoCount[23][$i];
					$indAutosomeMissCount = $indGenomeMissCount[$i] - $indChrMissCount[23][$i];
					$indAutosomeCountHomoAA = $indGenomeCountHomoAA[$i] - $indChrCountHomoAA[23][$i];
					$indAutosomeCountHomoBB = $indGenomeCountHomoBB[$i] - $indChrCountHomoBB[23][$i];
					$indAutosomeCountHeteroAB = $indGenomeCountHeteroAB[$i] - $indChrCountHeteroAB[23][$i];

					$indAutosomeExpHomoAA[$i] = $indPopGenomeSumHomoAA[$i] - $indPopChrSumHomoAA[23][$i];
					$indAutosomeExpHomoBB[$i] = $indPopGenomeSumHomoBB[$i] - $indPopChrSumHomoBB[23][$i];
					
					if ($ped[$i][4] == 1) {  #if male and chr23 then do different sum for hetero

						$indAutosomeExpHeteroAB[$i] = $indPopGenomeSumHeteroAB[$i];
					}
					else {

						$indAutosomeExpHeteroAB[$i] = $indPopGenomeSumHeteroAB[$i] - $indPopChrSumHeteroAB[23][$i];
					}

					$chiSqAutosome[$i] = (($indAutosomeCountHomoAA - $indAutosomeExpHomoAA[$i]) ** 2 / $indAutosomeExpHomoAA[$i]) + (($indAutosomeCountHeteroAB - $indAutosomeExpHeteroAB[$i]) ** 2 / $indAutosomeExpHeteroAB[$i]) + (($indAutosomeCountHomoBB - $indAutosomeExpHomoBB[$i]) ** 2 / $indAutosomeExpHomoBB[$i]);

					if ($popHetTest eq "verbose" || $popHetTest eq "perChrVerbose") {

						printf ("\t%7d\tObserve\t%7d\t%7d\t%7d\t", $indAutosomeMissCount, $indAutosomeCountHomoAA, $indAutosomeCountHeteroAB, $indAutosomeCountHomoBB);
						printf OUT ("\t%7d\tObserve\t%7d\t%7d\t%7d\t", $indAutosomeMissCount, $indAutosomeCountHomoAA, $indAutosomeCountHeteroAB, $indAutosomeCountHomoBB);
					}

					elsif ($popHetTest eq "perChr") {

						printf ("\t%7d\t%7d\t%7d\t%7d\t%7d\t%7.2f", $indAutosomeGenoCount, $indAutosomeMissCount, $indAutosomeCountHomoAA, $indAutosomeCountHeteroAB, $indAutosomeCountHomoBB, $chiSqAutosome[$i]);
						printf OUT ("\t%7d\t%7d\t%7d\t%7d\t%7d\t%7.2f", $indAutosomeGenoCount, $indAutosomeMissCount, $indAutosomeCountHomoAA, $indAutosomeCountHeteroAB, $indAutosomeCountHomoBB, $chiSqAutosome[$i]);
					}
				}
			}
		}

		if ($popHetTest eq "verbose" || $popHetTest eq "perChrVerbose") {
		
			print "\n\t\t\t\t\t";
			print OUT "\n\t\t\t\t\t";

			for ($i=0; $i<$pedno; $i++) {

				for($j=0; $j<=$#which; $j++){ # genotyped individual

					if($which[$j]==$whichsamples[$i]){

						printf ("\t\tExpect\t%7.1f\t%7.1f\t%7.1f\t%7.2f", $indAutosomeExpHomoAA[$i], $indAutosomeExpHeteroAB[$i], $indAutosomeExpHomoBB[$i], $chiSqAutosome[$i]);
						printf OUT ("\t\tExpect\t%7.1f\t%7.1f\t%7.1f\t%7.2f", $indAutosomeExpHomoAA[$i], $indAutosomeExpHeteroAB[$i], $indAutosomeExpHomoBB[$i], $chiSqAutosome[$i]);
					}
				}
			}					
		}

		print "\n";

		if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {
		
			print OUT "\n";
		}
	}

	$popGenomeMeanHet = $popGenomeSumHeteroAB / $genomeCountSNPs;

	if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

		print "-------\t-------\t-------\t-------\t-------\t-------";
		print OUT "-------\t-------\t-------\t-------\t-------\t-------";

		for ($i=0; $i<$pedno; $i++){

			for($j=0;$j<=$#which;$j++){ # genotyped individual

				if($which[$j]==$whichsamples[$i]){

					print "\t-------\t-------\t-------\t-------\t-------\t-------";
					print OUT "\t-------\t-------\t-------\t-------\t-------\t-------";
				}
			}
		}

		printf ("\ngenome\t%7d\t%7.4f\t%7.1f\t%7.1f\t%7.1f", $genomeCountSNPs, $popGenomeMeanHet, $popGenomeSumHomoAA, $popGenomeSumHeteroAB, $popGenomeSumHomoBB);
		printf OUT ("\ngenome\t%7d\t%7.4f\t%7.1f\t%7.1f\t%7.1f", $genomeCountSNPs, $popGenomeMeanHet, $popGenomeSumHomoAA, $popGenomeSumHeteroAB, $popGenomeSumHomoBB);
	}

	for ($i=0; $i<$pedno; $i++){

		for ($j=0; $j<=$#which; $j++){	

			if($which[$j] == $whichsamples[$i]){

				$chiSqGenome[$i] = (($indGenomeCountHomoAA[$i] - $indPopGenomeSumHomoAA[$i]) ** 2 / $indPopGenomeSumHomoAA[$i]) + (($indGenomeCountHeteroAB[$i] - $indPopGenomeSumHeteroAB[$i]) ** 2 / $indPopGenomeSumHeteroAB[$i]) + (($indGenomeCountHomoBB[$i] - $indPopGenomeSumHomoBB[$i]) ** 2 / $indPopGenomeSumHomoBB[$i]);

				if ($popHetTest eq "verbose" || $popHetTest eq "perChrVerbose") {

					printf ("\t%7d\tObserve\t%7d\t%7d\t%7d\t", $indGenomeMissCount[$i], $indGenomeCountHomoAA[$i], $indGenomeCountHeteroAB[$i], $indGenomeCountHomoBB[$i]);
					printf OUT ("\t%7d\tObserve\t%7d\t%7d\t%7d\t", $indGenomeMissCount[$i], $indGenomeCountHomoAA[$i], $indGenomeCountHeteroAB[$i], $indGenomeCountHomoBB[$i]);
				}

				elsif ($popHetTest eq "perChr") {

					printf ("\t%7d\t%7d\t%7d\t%7d\t%7d\t%7.2f", $indGenomeGenoCount[$i], $indGenomeMissCount[$i], $indGenomeCountHomoAA[$i], $indGenomeCountHeteroAB[$i], $indGenomeCountHomoBB[$i], $chiSqGenome[$i]);
					printf OUT ("\t%7d\t%7d\t%7d\t%7d\t%7d\t%7.2f", $indGenomeGenoCount[$i], $indGenomeMissCount[$i], $indGenomeCountHomoAA[$i], $indGenomeCountHeteroAB[$i], $indGenomeCountHomoBB[$i], $chiSqGenome[$i]);
				}
			}
		}
	}

	if ($popHetTest eq "verbose" || $popHetTest eq "perChrVerbose") {
		
		print "\n\t\t\t\t\t";
		print OUT "\n\t\t\t\t\t";

		for ($i=0; $i<$pedno; $i++) {

			for($j=0; $j<=$#which; $j++){ # genotyped individual

				if($which[$j]==$whichsamples[$i]){

					printf ("\t\tExpect\t%7.1f\t%7.1f\t%7.1f\t%7.2f", $indPopGenomeSumHomoAA[$i], $indPopGenomeSumHeteroAB[$i], $indPopGenomeSumHomoBB[$i], $chiSqGenome[$i]);
					printf OUT ("\t\tExpect\t%7.1f\t%7.1f\t%7.1f\t%7.2f", $indPopGenomeSumHomoAA[$i], $indPopGenomeSumHeteroAB[$i], $indPopGenomeSumHomoBB[$i], $chiSqGenome[$i]);
				}
			}
		}
	}

	if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

		print "\n-------\t-------\t-------\t-------\t-------\t-------";
		print OUT "\n-------\t-------\t-------\t-------\t-------\t-------";

		for ($i=0; $i<$pedno; $i++){

			for($j=0;$j<=$#which;$j++){ # genotyped individual

				if($which[$j]==$whichsamples[$i]){

					print "\t-------\t-------\t-------\t-------\t-------\t-------";
					print OUT "\t-------\t-------\t-------\t-------\t-------\t-------";
				}
			}
		}

		print "\nFor samples with missing genotypes, the expected counts are adjusted accordingly with respect to the allele frequencies for the specific missing SNP data.";
		print OUT "\nFor samples with missing genotypes, the expected counts are adjusted accordingly with respect to the allele frequencies for the specific missing SNP data.";

		if ($popHetTest eq "perChr") {

			print "\nUsing \"-popHetTest verbose\" or \"-popHetTest perChrVerbose\" will show sample specific expected counts.";
			print OUT "\nUsing \"-popHetTest verbose\" or \"-popHetTest perChrVerbose\" will show sample specific expected counts.";
		}

		if ($popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

			print "\nFor males, observed chr23 (X chromosome) heterozygous counts are shown only for reference, and not used in any calculation whatsoever.";
			print OUT "\nFor males, observed chr23 (X chromosome) heterozygous counts are shown for reference only, and not used in any calculation whatsoever.\n";
		}

		print "\nIndividuals with truncated names can be seen in full by opening the tab-delimited popHetTestDetails.txt table.\n";
	}

	open(OUT2, ">$tablesDir/popHetTest.txt") or print_usage "Cannot open file $tablesDir/popHetTest.txt, $!";

	print "\n\t\t\tchi^2 for the:";
	print "\nFamily\tMember\tSex\tGenome\tAutosome";
	print "\n-------\t-------\t-------\t-------\t--------";

	print OUT2 "Family\tMember\tSex\tGenomeChi^2\tAutosomeChi^2";

	for ($i=0; $i<$pedno; $i++){

		for($j=0;$j<=$#which;$j++){ # genotyped individual

			if($which[$j]==$whichsamples[$i]){

				printf ("\n%s\t%s\t%s\t%7.3f\t%8.3f", substr($ped[$i][0], 0, 7), substr($ped[$i][1], 0, 7), $sex[$i], $chiSqGenome[$i], $chiSqAutosome[$i]);
				printf OUT2 ("\n%s\t%s\t%s\t%7.3f\t%8.3f", $ped[$i][0], $ped[$i][1], $sex[$i], $chiSqGenome[$i], $chiSqAutosome[$i]);
				
				if ( length($ped[$i][0]) > 7) {

					print "\tFamily name truncated.";
				}
				
				if ( length($ped[$i][0]) > 7) {
				
					print "\tMember's name truncated.";
				}
				
				if ( length($ped[$i][0]) > 7 || length($ped[$i][1]) > 7) {

					print "\tSee the tab-delimited popHetTest.txt table for full family/member names.";
				}
			}
		}
	}

	print "\n-------\t-------\t-------\t-------\t--------\n";

	if (defined($randomSNP)) {
	
		print "\nWARNING: Random SNPs selected from each bin.  This will adversely effect the above analysis.  Re-run LINKDATAGEN without the \"-randomSNP\" option.\n";
		print OUT2 "\nWARNING: Random SNPs selected from each bin.  This will adversely effect the above analysis.  Re-run LINKDATAGEN without the \"-randomSNP\" option.\n";

		if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {
		
			print OUT "\n\nWARNING: Random SNPs selected from each bin.  This will adversely effect the above analysis.  Re-run LINKDATAGEN without the \"-randomSNP\" option.\n";
		}
	}

	if (defined($removeWFHBS)) {
	
		print "\nWARNING: \"Within-family homozygosity-by-state\" SNP markers have been removed.  This will make the available SNPs more heterozygous than by chance alone, thus adversely effecting this analysis.  It is recommended NOT to use the \"-removeWFHBS\" option when running this analysis.\n";
		print OUT2 "\nWARNING: \"Within-family homozygosity-by-state\" SNP markers have been removed.  This will make the available SNPs more heterozygous than by chance alone, thus adversely effecting this analysis.  It is recommended NOT to use the \"-removeWFHBS\" option when running this analysis.\n";

		if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

			print OUT "\n\nWARNING: \"Within-family homozygosity-by-state\" SNP markers have been removed.  This will make the available SNPs more heterozygous than by chance alone, thus adversely effecting this analysis.  It is recommended NOT to use the \"-removeWFHBS\" option when running this analysis.\n";
		}
	}

	if ($popHetTest eq "verbose" || $popHetTest eq "perChr" || $popHetTest eq "perChrVerbose") {

		close(OUT);
	}

	close(OUT2);
	
	return(\@chiSqAutosome);
}

sub print_genos_test{
my $i;
my $k;
my $keycnt=0;
my @testsnp=();
	for($i=0;$i<=2;$i++){
		$testsnp[$i]=$chr_order_snp[1][$i][0];
	}
	for($i=0;$i<=$#testsnp;$i++){
		print "TEST SNP $testsnp[$i]\n";
		print "genos\n";	
		for($k=0;$k<=$#which;$k++){
			print "$genos{$testsnp[$i]}[$k] ";
		}	
		print "\n";
		print "genos_orig\n";	
		for($k=0;$k<=$#which;$k++){
			print "$genos_orig{$testsnp[$i]}[$k] ";
		}	
		print "\n";
	}	
}

# ----------------OUTPUT functions ------------------------------------------
# merlin output functions for .ped, .freq, .dat and .map files (QTDT format)
sub print_pre_merlin{
my $i;
my $j;
my $k;
my $l;
my $filename;
my $geno;
my $c;
my $chrom;

	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$merlinDir."/merlin_".$chrom."_".$pfile.".ped";
		#print "writing merlin .ped file $filename...\n";
		open(OUT,">$filename") or print_usage "Cannot open file $filename, $!";
		IND: for($i=0;$i<$pedno;$i++){
			for($j=0;$j<=$#which;$j++){ # genotyped individual
				if($which[$j]==$whichsamples[$i]){ # Why do I need this line at all?
					for($k=0;$k<5;$k++){

						print OUT "$ped[$i][$k]\t";
					}
					#print "Current chromo is ",$c+1,". Snp count is $chr_snp_cnt[$c]\n";
					for($k=0;$k<$chr_snp_cnt[$c];$k++){

						$geno=$genos{$chr_order_snp[$c][$k][0]}[$j];

						if (!defined($actg)) {

							if($geno== -1){
								print OUT "0/0\t";
							}
							if($geno==0){
								print OUT "1/1\t";
							}
							if($geno==1){
								print OUT "1/2\t";
							}
							if($geno==2){
								print OUT "2/2\t";
							}
						}
						else {

							if($geno== -1){
								print OUT "0/0\t";
							}
							if($geno==0){
								print OUT "1/1\t";
							}
							if($geno==1){
								print OUT "1/2\t";
							}
							if($geno==2){
								print OUT "2/2\t";
							}
						}						
						if($geno==3){
							print_usage("Problem with b/crlmm file. Assumed to be brlmm but appears to be crlmm. Check, and switch off -crlmm option if necessary.\n");
						}
						if(!defined($geno)){
							print_usage("Undefined genotype in merlin subroutine\nchr = $c, ind $i, SNP $k\n");
						}
					}
					print OUT "$ped[$i][5]\n";
					next IND;
				}
			}
			# not genotyped, all missing..
			for($k=0;$k<5;$k++){
				print OUT "$ped[$i][$k]\t";
			}
			for($k=0;$k<$chr_snp_cnt[$c];$k++){
				print OUT "0/0\t";
			}
			print OUT "$ped[$i][5]\n";
		}

		close(OUT);
	}

	open (OUTG, ">$merlinDir"."/genome/merlin_genome_".$pfile.".ped") || die print ("Can't open file ".$merlinDir."/genome/merlin_genome_".$pfile.".ped");

	IND: for($i=0;$i<$pedno;$i++){

		for($k=0;$k<6;$k++){

			print OUTG "$ped[$i][$k]\t";
		}

		for($j=0;$j<=$#which;$j++){ # genotyped individual

			if($which[$j]==$whichsamples[$i]){

				for($c=0;$c<$numchrom;$c++){

					if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

						if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {

							next;
						}
					}

					$chrom=$c+1;

					for($k=0;$k<$chr_snp_cnt[$c];$k++){

						$geno=$genos{$chr_order_snp[$c][$k][0]}[$j];

						if (!defined($actg)) {

							if($geno== -1){
								print OUTG "0/0\t";
							}
							if($geno==0){
								print OUTG "1/1\t";
							}
							if($geno==1){
								print OUTG "1/2\t";
							}
							if($geno==2){
								print OUTG "2/2\t";
							}
						}
						else {

							if($geno== -1){
								print OUTG "0/0\t";
							}
							if($geno==0){
								print OUTG "1/1\t";
							}
							if($geno==1){
								print OUTG "1/2\t";
							}
							if($geno==2){
								print OUTG "2/2\t";
							}
						}
						if($geno==3){
							print_usage("Problem with b/crlmm file. Assumed to be brlmm but appears to be crlmm. Check, and switch off -crlmm option if necessary.\n");
						}
						if(!defined($geno)){
							print_usage("Undefined genotype in merlin subroutine\nchr = $c, ind $i, SNP $k\n");
						}
					}
				}

				print OUTG "\n";
				next IND;
			}
		}

		# not genotyped, all missing

		for($c=0;$c<$numchrom;$c++) {

			if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ) {    #test if regions are being used and if so whether region included

				if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {

					next;
				}
			}

			$chrom=$c+1;

			for($k=0;$k<$chr_snp_cnt[$c];$k++){

				print OUTG "0/0\t";
			}
		}

		print OUTG "\n";
	}

	close(OUTG);
}

sub print_dat_merlin{
my $c;
my $j;
my $chrom;
my $filename;

	open (OUTA, ">$merlinDir"."/genome/merlin_autosome_".$pfile.".dat") || die print ("Can't open file ".$merlinDir."/genome/merlin_autosome_".$pfile.".dat");
	open (OUTX, ">$merlinDir"."/genome/merlin_X_".$pfile.".dat") || die print ("Can't open file ".$merlinDir."/genome/merlin_X_".$pfile.".dat");
	print OUTA "A Disease\n";
	print OUTX "A Disease\n";

	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$merlinDir."/merlin_".$chrom."_".$pfile.".dat";
		#print "writing merlin .dat file $filename...\n";
		open(OUT,">$filename")|| die print "Can't open file $filename\n";
		
		if ($chrom < 23) {

			for($j=0;$j<$chr_snp_cnt[$c];$j++){

				print OUT "M $chr_order_snp[$c][$j][0]\n";
				print OUTA "M $chr_order_snp[$c][$j][0]\n";
				print OUTX "S2 $chr_order_snp[$c][$j][0]\n";
			}
		}

		if ($chrom == 23) {

			for($j=0;$j<$chr_snp_cnt[$c];$j++){

				print OUT "M $chr_order_snp[$c][$j][0]\n";
				print OUTA "S2 $chr_order_snp[$c][$j][0]\n";
				print OUTX "M $chr_order_snp[$c][$j][0]\n";
			}
		}
		
		print OUT "A disease\n";
		close(OUT);
	}
	
	close OUTA;
	close OUTX;
}

sub print_map_merlin{
my $c;
my $j;
my $chrom;
my $sameGeneticPos;
my $filename;

	open (OUTG, ">$merlinDir"."/genome/merlin_genome_".$pfile.".map") || die print ("Can't open file ".$merlinDir."/merlin_genome_".$pfile.".map");
	print OUTG "CHROMOSOME\tMARKER\tPOSITION\n";

	for($c=0;$c<$numchrom;$c++){

		$sameGeneticPos=0;

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$merlinDir."/merlin_".$chrom."_".$pfile.".map";

		#print "writing merlin .map file $filename...\n";

		open(OUT,">$filename") || die print "Can't open file $filename\n";
		print OUT "CHROMOSOME\tMARKER\tPOSITION\n";

		for($j=0;$j<$chr_snp_cnt[$c];$j++){

			if ($chr_order_snp[$c][$j][3] != $chr_order_snp[$c][$j-1][3]) {

				$sameGeneticPos=0;

				printf OUT ("%s\t%s\t%1.4f\n", $chr_order_snp[$c][$j][1], $chr_order_snp[$c][$j][0], $chr_order_snp[$c][$j][3]);
				printf OUTG ("%s\t%s\t%1.4f\n", $chr_order_snp[$c][$j][1], $chr_order_snp[$c][$j][0], $chr_order_snp[$c][$j][3]);

			}
			else {

				printf OUT ("%s\t%s\t%1.4f%04d\n", $chr_order_snp[$c][$j][1], $chr_order_snp[$c][$j][0], $chr_order_snp[$c][$j][3], ++$sameGeneticPos);
				printf OUTG ("%s\t%s\t%1.4f%04d\n", $chr_order_snp[$c][$j][1], $chr_order_snp[$c][$j][0], $chr_order_snp[$c][$j][3], $sameGeneticPos);
			}
		}
		close OUT;
	}
	
	close OUTG;
}

# prints out all merlin freq files, one per chromo.
sub print_freq_merlin{
my $c;
my $j;
my $snpfreq;
my $chrom;
my $filename;

	open (OUTG, ">$merlinDir"."/genome/merlin_genome_".$pfile.".freq") || die print ("Can't open file ".$merlinDir."/merlin_genome_".$pfile.".freq");

	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$merlinDir."/merlin_".$chrom."_".$pfile.".freq";
		#print "writing merlin .freq file $filename...\n";
		open(OUT,">$filename") || die print "Can't open file $filename\n";
		for($j=0;$j<$chr_snp_cnt[$c];$j++){

			print OUT "M $chr_order_snp[$c][$j][0]\n";
			print OUTG "M $chr_order_snp[$c][$j][0]\n";

			$snpfreq=$chr_order_snp[$c][$j][4];
			#if($snpfreq>0 && $snpfreq<1){

				printf OUT ("F %1.3f %1.3f\n", $snpfreq+0.001, 1-$snpfreq+0.001);
				printf OUTG ("F %1.3f %1.3f\n", $snpfreq+0.001, 1-$snpfreq+0.001);
		}
		close(OUT);
	}

	close OUTG;
}

# making the 23 merlin chromosome run files.
sub print_merlin_run_files{
my $i;
my $filename;
my $whichmerlin;
	$whichmerlin="merlin";
	for($i=1;$i<=$numchrom;$i++){

		if ( (%regions && !defined($regions{$i})) || !defined($chr_snp_cnt_orig[$i-1]) ){    #test if regions are being used and if so whether region included

			if ($i <= 22 || ($i == 23 && !defined($regions{"X"}))) {	###note uses i not c and is one extra

				next;
			}
		}


		$filename=$merlinDir."/merlin_".$i."_".$pfile.".in";
		open(OUT,">$filename") || die print "Can't open file $filename\n";
		if($i==23){
			$whichmerlin="minx";
		}
		# all necessary input files.
		# make a special case for the X chromosome - not yet done.
		# getting all double recombinant errors first.
		print OUT "mv param.tbl temp.tbl\n"; # so parametric linkage analysis isn't run in the error checking (double recombinant check) step that is next
		print OUT "$whichmerlin -d merlin_",$i,"_",$pfile,".dat -p merlin_",$i,"_",$pfile,".ped -f merlin_",$i,"_",$pfile,".freq -m merlin_",$i,"_",$pfile,".map --prefix merlin_",$i,"_",$pfile;
		# analysis options wanted
		print OUT " --error\n";
		# print OUT "cp merlin.err merlin_",$i,"_",$pfile,".err\n";  # not required as now using --prefix
		# running pedwipe to get rid of double recombinant errors
		print OUT "pedwipe -d merlin_",$i,"_",$pfile,".dat -p merlin_",$i,"_",$pfile,".ped -e merlin_",$i,"_",$pfile,".err\n";
		print OUT "mv wiped.dat merlin_",$i,"_",$pfile,".dat\n";
		print OUT "mv wiped.ped merlin_",$i,"_",$pfile,".ped\n";
		print OUT "mv temp.tbl param.tbl\n";
		print OUT "$whichmerlin -d merlin_",$i,"_",$pfile,".dat -p merlin_",$i,"_",$pfile,".ped -f merlin_",$i,"_",$pfile,".freq -m merlin_",$i,"_",$pfile,".map";
		# analysis options wanted
		#print OUT " --ibd --extended --founders --best --error --model param.tbl --pdf > merlin_",$i,"_",$pfile,".out";
		# add in --npl for all calculations
		print OUT " --smallswap --megabytes:9999 --founders --grid 0.3 --pairs --ibd --extended --best --exp --model param.tbl --pdf --tabulate --prefix merlin_",$i,"_",$pfile,"> merlin_",$i,"_",$pfile,".out\n";
		print OUT "rm wiped.freq\n";
		close(OUT);
	}
	
		open (OUTG, ">$merlinDir"."/genome/merlin_genome_".$pfile.".in") || die print ("Can't open file ".$merlinDir."/merlin_genome_".$pfile.".in");

		print OUTG "mv param.tbl temp.tbl\n"; # so parametric linkage analysis isn't run in the error checking (double recombinant check) step that is next
		print OUTG "merlin -d merlin_autosome_",$pfile,".dat -p merlin_genome_",$pfile,".ped -f merlin_genome_",$pfile,".freq -m merlin_genome_",$pfile,".map --prefix merlin_autosome_",$pfile," --error\n";
		print OUTG "pedwipe -d merlin_autosome_",$pfile,".dat -p merlin_genome_",$pfile,".ped -e merlin_autosome_",$pfile,".err\n";
		print OUTG "mv wiped.dat merlin_autosome_wiped_",$pfile,".dat\n";
		print OUTG "mv wiped.ped merlin_autosome_wiped_",$pfile,".ped\n";
		print OUTG "minx -d merlin_X_",$pfile,".dat -p merlin_genome_",$pfile,".ped -f merlin_genome_",$pfile,".freq -m merlin_genome_",$pfile,".map --prefix merlin_X_",$pfile," --error\n";
		print OUTG "pedwipe -d merlin_X_",$pfile,".dat -p merlin_genome_",$pfile,".ped -e merlin_X_",$pfile,".err\n";
		print OUTG "mv wiped.dat merlin_X_wiped_",$pfile,".dat\n";
		print OUTG "mv wiped.ped merlin_X_wiped_",$pfile,".ped\n";
		print OUTG "rm wiped.freq\n";
		print OUTG "mv temp.tbl param.tbl\n";
		print OUTG "merlin -d merlin_autosome_wiped_",$pfile,".dat -p merlin_autosome_wiped_",$pfile,".ped -f merlin_genome_",$pfile,".freq -m merlin_genome_",$pfile,".map";
		# analysis options wanted
		print OUTG " --smallswap --megabytes:9999 --founders --grid 0.3 --pairs --ibd --extended --best --exp --model param.tbl --pdf --tabulate --prefix merlin_autosome_",$pfile," > merlin_autosome_",$pfile,".out\n";
		print OUTG "minx -d merlin_X_wiped_",$pfile,".dat -p merlin_X_wiped_",$pfile,".ped -f merlin_genome_",$pfile,".freq -m merlin_genome_",$pfile,".map";
		# analysis options wanted
		print OUTG " --smallswap --megabytes:9999 --founders --grid 0.3 --pairs --ibd --extended --best --exp --model param.tbl --pdf --tabulate --prefix merlin_X_",$pfile," > merlin_X_",$pfile,".out\n";
		close OUTG;

		open (OUTR, ">$merlinDir"."/README.txt") || die print ("Can't open file ".$merlinDir."/README.txt");
		print OUTR "This folder contains most of the files required to run MERLIN.\n";
		print OUTR "Please refer to the MERLIN tutorial if you are not familiar with using MERLIN or linkage analysis:\n";
		print OUTR "http://www.sph.umich.edu/csg/abecasis/merlin/tour/linkage.html\n\n";
		print OUTR "***  Before continuing, you should make a param.tbl file \(see below\).  ***\n\n";
		print OUTR "The files in this folder include merlin_*.ped|dat|freq|map|in files.\n";
		print OUTR "There is a set of each files for each chromosome \(1 to 22, and X\).\n";
		print OUTR "You can use these files to perform linkage analysis on individual chromosomes.\n";
		print OUTR "There is also a set of merlin_genome_*.ped|dat|freq|map|in files in the \"genome\" folder.\n\n";
		print OUTR "The *.in files are scripts that contain all the commands required to run MERLIN.\n";
		print OUTR "If you are using a UNIX-based system, you might try typing \"bash merlin_*.in\" filling in * as necessary.\n";
		print OUTR "Alternatively, you can try copying and pasting the lines of the *.in files to run MERLIN.\n\n";
		print OUTR "To use our *.in scripts, it is necessary for you to create a file called \"param.tbl\" that is required by MERLIN.\n";
		print OUTR "This file describes the disease model (dominant, recessive etc...) that you are testing for linkage in your dataset.\n";
		print OUTR "In the simplest form, the \"param.tbl\" file consists of a single line. Two examples are given below that you may use:\n\n";
		print OUTR "Disease	0.001	0.0001,1.0,1.0	Rare_dominant\n\n";
		print OUTR "Disease	0.001	0.0001,0.0001,1.0	Rare_recessive\n\n";
		print OUTR "The first value \"Disease\" is required to match up with the  \"Affection Status\" column in the *.ped files, as indicated by the *.dat files.\n";
		print OUTR "The second value, \"0.001\", is the disease allele frequency in the general population.\n";
		print OUTR "The next set of values are the penetrances of the disease for carrying 0, 1 or 2 copies of the disease allele.\n";
		print OUTR "That is, the penetrances of the disease for being homozygous \(wildtype\), heterozygous, or homozygous \(disease allele\).\n";
		print OUTR "The final value is simply a name for your disease model.\n\n";
		print OUTR "***  If you wish to try different models \(dominant, recessive etc...\) then do NOT repeatedly run the *.in files.  ***\n\n";
		print OUTR "Doing so will re-run the error detection option \(--error\) of MERLIN which is not desirable.\n\n";
		print OUTR "Beginners should create different folders of *.ped|dat|freq|map|in files generated by LINKDATAGEN for each disease model tested.\n";
		print OUTR "Experienced/advanced users should run MERLIN once with the --error option, and then subsequently without this option.\n\n";
		print OUTR "Alternatively (and better still), the \"param.tbl\" file need not be a single line; alternative models can be included on separate lines.\n\n";
		print OUTR "If at all unsure, then re-run LINKDATAGEN to create a new folder of MERLIN input files each time you wish to test a different model.\n";
		close OUTR;

		open (OUTR, ">$merlinDir"."/genome/README.txt") || die print ("Can't open file ".$merlinDir."/genome/README.txt");
		print OUTR "Remember to create a \"param.tbl\" file prior to running the merlin_genome_*.in script.\n";
		print OUTR "In the simplest form it consists of a single line.  Two examples are given below that you may use:\n\n";
		print OUTR "Disease	0.001	0.0001,1.0,1.0	Rare_dominant\n\n";
		print OUTR "Disease	0.001	0.0001,0.0001,1.0	Rare_recessive\n\n";
		print OUTR "The first value \"Disease\" is required to match up with the \"Affection Status\" column in the *.ped files, as indicated by the *.dat files.\n";
		print OUTR "The second value, \"0.001\", is the disease allele frequency in the general population.\n";
		print OUTR "The next set of values are the penetrances of the disease for carrying 0, 1 or 2 copies of the disease allele.\n";
		print OUTR "That is, the penetrances of the disease for being homozygous \(wildtype\), heterozygous, or homozygous \(disease allele\).\n";
		print OUTR "The final value is simply a name for your disease model.\n\n";
		print OUTR "There are two different *.dat files in this folder.  This is because the X chromosome data is analysed separately with MINX, rather than MERLIN.\n";
		print OUTR "Once you have created a \"param.tbl\" file you should be able to run the merlin_genome_*.in script.\n\n";
		print OUTR "***  If you wish to try different models \(dominant, recessive etc...\) then do NOT repeatedly run the *.in file.  ***\n\n";
		print OUTR "Doing so will re-run the error detection option \(--error\) of MERLIN which is not desirable.\n";
		print OUTR "See the README.txt file in the above folder.\n";
		close OUTR;
}

# Allegro output files: 
# pre file for allegro, matched LINKAGE format
sub print_pre_allegro{
my $i;
my $j;
my $k;
my $filename;
my $geno;
my $c;
my $chrom;
my $pedidcnt=0;
	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$allegroDir."/allegro_".$chrom."_".$pfile.".pre";
		#print "writing allegro .pre file $filename...\n";
		open(OUT,">$filename") or print_usage "Cannot open file $filename, $!";
		IND: for($i=0;$i<$pedno;$i++){
			$pedidcnt=0;
			for($j=0;$j<=$#whichsamples;$j++){ # genotyped individual
				if($i==$j && $whichsamples[$i]!=0){
					#print "individual $ped[$i][1] in pedigree, corresponds to $whichsamples[$j]\n";
					#print "pedidcnt=$pedidcnt\n";
					for($k=0;$k<=5;$k++){
						print OUT "$ped[$i][$k] ";
					}
					for($k=0;$k<$chr_snp_cnt[$c];$k++){
						#print "Current chromo is ",$c+1,". Snp count is $chr_snp_cnt[$c]\n";
						#print "SNP is $chr_order_snp[$c][$k][0], ind count in geno array is $pedidcnt\n";
						$geno=$genos{$chr_order_snp[$c][$k][0]}[$pedidcnt];
						if($geno== -1){ 
							print OUT "0 0\t";
						}
						if($geno==0){ 
							print OUT "1 1\t";
						}
						if($geno==1){ 
							print OUT "1 2\t";
						}
						if($geno==2){
							print OUT "2 2\t";
						}
						if($geno==3){
							print_usage("Problem with BRLMM calls. Found an impossible \"3\" genotype call. Check that the -crlmm option is on.");
						}
					}
					print OUT "\n";
					next IND;
				}
				if($whichsamples[$j]!=0){ # statement A
					$pedidcnt+=1;
				}
			}
			# not genotyped, all missing..
			for($k=0;$k<=5;$k++){
				print OUT "$ped[$i][$k] ";
			}
			for($k=0;$k<$chr_snp_cnt[$c];$k++){
				print OUT "0 0 ";
			}
			print OUT "\n";
		}
		close(OUT);
	}
}

# allegro style .dat files which are the LINKAGE style .dat files.
sub print_dat_allegro{
my $j;
my $c;
my $chrom;
my $filename;
my $snpfreq;
my $dist;
	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$allegroDir."/allegro_".$chrom."_".$pfile.".dat";
		#print "writing allegro .dat file $filename...\n";
		open(OUT,">$filename") or print_usage "Cannot write file $filename, $!";
		print OUT $chr_snp_cnt[$c]+1," 0 0 5\n"; 
		print OUT "0 0.0 0.0 0\n"; 
		for($j=1;$j<=$chr_snp_cnt[$c];$j++){
			print OUT "$j ";
		}
		print OUT "\n";
		print OUT "1 2\n";
		# disease alle freqs
		print OUT "0.9999 0.0001\n"; # CHANGED BY USER AS NECESSARY
		print OUT "1\n"; 
		print OUT "0.0001 0.9999 0.9999\n";# CHANGED BY USER AS NECESSARY
		# now all the marker data:
		for ($j=0;$j<$chr_snp_cnt[$c];$j++){
			print OUT "3 2 # $chr_order_snp[$c][$j][0]\n"; # contains marker name as well
			$snpfreq=$chr_order_snp[$c][$j][4];
			if($snpfreq>0 && $snpfreq<1){
				printf OUT "%1.2f",$snpfreq;
				print OUT " ";
				printf OUT "%1.2f",1-$snpfreq;
				print OUT "\n";
			}
			else{ # zero freq in hapmap sample
				if($snpfreq==0){
					print OUT "0.01 0.99\n";
				}
				else{
					print OUT "0.99 0.01\n";
				}
			}
		}
		# now the map stuff
		print OUT "0 0\n";
		for($j=1;$j<$chr_snp_cnt[$c];$j++){
			$dist=$chr_order_snp[$c][$j][3]-$chr_order_snp[$c][$j-1][3];
			printf OUT "%1.2f",$dist;
			print OUT " ";
		}
		print OUT "\n";
		print OUT "1 0.10000 0.450000\n";
		close(OUT);
	}
}

# making the 23 allegro chromosome run files.
# unlike the merlin files above these still have to be run with the command allegro allegro_x.in in a separate
# shell file 
sub print_allegro_run_files{
my $i;
my $haplofile;
my $ihaplofile;
my $founderfile;
my $inherfile;
my $paramfile;
my $uninffile;
my $filename;
	$filename=$allegroDir."/runallegrofiles_".$pfile.".sh";
	open(OUT1,">$filename") or print_usage "Cannot write file $filename, $!";
	for($i=1;$i<=$numchrom;$i++){

		if (%regions && !defined($regions{$i})){    #test if regions are being used and if so whether region included

			if ($i <= 22 || ($i == 23 && !defined($regions{"X"}))) {	###note uses i not c and is one extra

				next;
			}
		}

		$filename=$allegroDir."/allegro_".$i."_".$pfile.".in";
		open(OUT,">$filename")  or print_usage "Cannot write file $filename, $!";
		print OUT1 "allegro allegro_".$i."_".$pfile.".in\n";
		# all necessary input files.
		# make a special case for the X chromosome
		print OUT "PREFILE allegro_".$i."_".$pfile.".pre\n"; 
		print OUT "DATFILE allegro_".$i."_".$pfile.".dat\n\n"; 
		$paramfile="param_".$i.".out";
		if($i<23){
			print OUT "MODEL mpt par $paramfile\n"; # can stick in model params here if desired
		}
		else{
			print OUT "MODEL mpt par X $paramfile\n"; # can stick in model params here if desired		
		}
		$haplofile="haplo_".$i."_".$pfile.".out";
		$ihaplofile="ihaplo_".$i."_".$pfile.".out";
		$founderfile="founder_".$i."_".$pfile.".out";
		$inherfile="inher_".$i."_".$pfile.".out";
		print OUT "HAPLOTYPE $haplofile $ihaplofile $founderfile $inherfile\n";
		$uninffile="uninf_".$i."_".$pfile.".out";	
		print OUT "UNINFORMATIVE $uninffile\n";
		print OUT "PAIRWISEIBD mpt genotyped\n"; # generates two files: prior.mpt and posterior.mpt. 	
		close(OUT);
	}	
	close(OUT1);
}

# PREST mcpeek & co, to check for sample swaps/pedigree errrors
# NOte that the pedigree needs to have all integer pedids and indids. I'm assuming for the moment that this is ok. More easily fixed ad hoc outside of linkdatagen for the moment.
# idx file is very similar to .dat file but has less in it. Needs recomb fractions, not cM distances.
sub print_PREST_idx{
my $j;
my $c;
my $chrom;
my $filename;
my $snpfreq;
my $dist;
my $newdist;
my $first;
my $second;
	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$prestDir."/prest_".$chrom."_".$pfile.".idx";

		#	print "writing PREST .idx file $filename...\n";

		open(OUT,">$filename")  or print_usage "Cannot write file $filename, $!";
		print OUT "$chr_snp_cnt[$c]\n"; 
		# now all the marker data:
		for ($j=0;$j<$chr_snp_cnt[$c];$j++){
			print OUT "3 2 # $chr_order_snp[$c][$j][6]\n"; # contains marker name as well
			$snpfreq=$chr_order_snp[$c][$j][4];
			$first=sprintf("%1.2f",$snpfreq);
			$second=1-$first;
			if($snpfreq>0 && $snpfreq<1){
				printf OUT "$first $second\n";
			}
			else{ # zero freq in hapmap sample
				if($snpfreq==0){
					print OUT "0.01 0.99\n";
				}
				else{
					print OUT "0.99 0.01\n";
				}
			}
		}
		# now the map stuff
		for($j=1;$j<$chr_snp_cnt[$c];$j++){
			$dist=$chr_order_snp[$c][$j][3]-$chr_order_snp[$c][$j-1][3];
			# change from cM to r via Haldane (lazy!) map function
			$newdist=(1-exp(-2*$dist/100))/2;
			printf OUT "%1.4f",$newdist;
			print OUT " ";
		}
		close(OUT);
	}
}

# output files for MORGAN (Thomson et al)
sub print_MORGAN{
my $c;
my $filename;
my $mfilename;
my $pfilename;
my $chrom;
my $i;
my $snpfreq;
my $first;
my $second;
my $j;
my $k;
my $geno;
my $cnt=0;
my $lmauto=1; # lmauto=0 if lm_auto otherwise lm_markers or lm_bayes is required.
	# print out the parameter file
	$filename=$morganDir."/runmorgan.sh";
	open(OUT3,">$filename") or print_usage "Cannot write file $filename, $!";
	open(OUT4,">$morganDir/input.seeds") or print_usage "Cannot write file $morganDir/input.seeds, $!";
	print OUT4 "set sampler seeds  0x46777249 0xcf3dbd61\n";
	print OUT4 "set sampler seeds  0xc28e1359 0x81239ee9\n";
	print OUT4 "set sampler seeds  0xc1cda0f9 0xc81f30ae\n";
	close(OUT4);
	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		print OUT3 "lm_markers morgan_".$chrom."_".$pfile.".par\n";
		$filename=$morganDir."/morgan_".$chrom."_".$pfile.".par";
		$mfilename=$morganDir."/morgan_".$chrom."_".$pfile.".mar";

		#	print "writing MORGAN .mar file $filename...\n";
		#	print "writing MORGAN .par file $mfilename...\n";

		$pfilename=$morganDir."/morgan_".$chrom."_".$pfile.".ped";

		#	print "writing MORGAN .ped file $pfilename...\n";

		open(OUT1,">$filename") or print_usage "Cannot write file $filename, $!";
		print OUT1 "input pedigree file 'morgan_".$chrom."_".$pfile.".ped'\n";
		print OUT1 "input seed file 'input.seeds'\n";
		print OUT1 "input pedigree size $pedno\n";
		print OUT1 "input pedigree record names 3 integers 2\n"; # two integers to describe a) gender (assumed to be present) and b) affectedness =  a binary trait
		print OUT1 "input marker data file 'morgan_".$chrom."_".$pfile.".mar'\n"; # specifies marker data file name
		# for lm_markers or lm_bayes
		if($lmauto==1){
			print OUT1 "map trait 1 all interval proportions 0.5\n";
		}
		# for lm_auto
		if($lmauto==0){
			print OUT1 "map trait ";
			for($i=0;$i<$chr_snp_cnt[$c];$i++){
				print OUT1 "1 ";
			}
			print OUT1 "marker ";
			for($i=1;$i<=$chr_snp_cnt[$c];$i++){
				print OUT1 "$i ";
			}
			print OUT1 "distances ";
			for($i=0;$i<$chr_snp_cnt[$c];$i++){
				printf OUT1 "%1.1f", $chr_order_snp[$c][$i][3];
				print OUT1 " ";
			}
			print OUT1 "\n";
		}
		print OUT1 "set trait data discrete\n";
		print OUT1 "set incomplete penetrances 0.05 0.9 0.9\n";
		print OUT1 "select all markers traits 1\n";
		print OUT1 "select all markers traits 1\n";
		print OUT1 "set trait 1 freqs 0.9999  0.0001\n";
		# Monte Carlo setup and requests
		print OUT1 "sample by scan\n";
		print OUT1 "set L-sampler probability 0.2\n";
		## For real analyses, recommended number of iterations is of order 10^4
		##   at each simulation position
		print OUT1 "set MC iterations 30000\n";
		print OUT1 "set burn-in iterations 1000\n";
		close(OUT1);
		# print out the pedigree file - no genotyping data allowed in this, and no pedigree qualifier in front of person's id
		open(OUT,">$pfilename") or print_usage "Cannot write file $pfilename, $!";
		print OUT "****************\n";
		for($i=0;$i<$pedno;$i++){
			for($k=1;$k<=5;$k++){
				print OUT "$ped[$i][$k] ";
			}
			print OUT "\n"; 
		}
		# print out the marker file containing the genotype data
		close(OUT);
		open(OUT2,">$mfilename") or print_usage "Cannot write file $mfilename, $!";
		# map data first
		#map marker positions 0.94 3.29 3.56 4.26 4.59 ... 275.64 275.69
		print OUT2 "map marker positions ";
		for($i=0;$i<$chr_snp_cnt[$c];$i++){
			printf OUT2 "%1.2f", $chr_order_snp[$c][$i][3];
			print OUT2 " ";
		}
		print OUT2 "\n\n";
		# marker frequency data
		#set markers 1 	freqs 0.85 0.15
		for ($j=0;$j<$chr_snp_cnt[$c];$j++){
			print OUT2 "set markers ",$j+1," freq ";
			$snpfreq=$chr_order_snp[$c][$j][4];
			$first=sprintf("%1.2f",$snpfreq);
			$second=1-$first;
			if($snpfreq>0 && $snpfreq<1){
				printf OUT2 "$first $second\n";
			}
			else{ # zero freq in hapmap sample
				if($snpfreq==0){
					print OUT2 "0.01 0.99\n";
				}
				else{
					print OUT2 "0.99 0.01\n";
				}
			}
		}
		print OUT2 "\n\n";		
		# genotyping data
		# format is indid m1geno m2geno etc
		print OUT2 "set markers $chr_snp_cnt[$c] data \n";
		$cnt=0;
		for($i=0;$i<$pedno;$i++){
			if($whichsamples[$i]!=0){
				print OUT2 "$ped[$i][1] ";
				for($k=0;$k<$chr_snp_cnt[$c];$k++){
					#print "Current chromo is ",$c+1,". Snp count is $chr_snp_cnt[$c]\n";
					#print "SNP is $chr_order_snp[$c][$k][0], ind count in geno array is $pedidcnt\n";
					$geno=$genos{$chr_order_snp[$c][$k][0]}[$cnt];
					if($geno== -1){ 
						print OUT2 "0 0\t";
					}
					if($geno==0){ 
						print OUT2 "1 1\t";
					}
					if($geno==1){ 
						print OUT2 "1 2\t";
					}
					if($geno==2){
						print OUT2 "2 2\t";
					}
				}
				$cnt+=1;
				print OUT2 "\n";
			}
		}
		close(OUT2);
	}
	close(OUT3);
}

# output files for PLINK
sub print_PLINK{
my $c;
my $chrom;
my $i;
my $j;
my $k;
my $geno;
my $temp;
my $cnt=0;
my $filename;
	$filename=$plinkDir."/plink.map";
	open(OUT2,">$filename") or print_usage "Cannot write file $filename, $!";
	# print out the parameter file
	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		for($i=0;$i<$chr_snp_cnt[$c];$i++){
			print OUT2 "$chrom ";
			print OUT2 "$chr_order_snp[$c][$i][6] ";
			printf OUT2 "%1.2f", $chr_order_snp[$c][$i][3];
			print OUT2 " ";
			#$temp=1000000*$chr_order_snp[$c][$i][3];
			printf OUT2 "%1.0f", $chr_order_snp[$c][$i][7];
			print OUT2 " ";
			print OUT2 "\n";
		}
	}
	close(OUT2);
	$filename=$plinkDir."/plink.ped";	
	open(OUT2,">$filename") or print_usage "Cannot write file $filename, $!";
	for($i=0;$i<$pedno;$i++){
		if( $whichsamples[$i]!=0 || ( defined($plOutput) && $plOutput eq "allSamples" ) ) {

			if ( defined($plOutput) && $plOutput eq "makeAllFounders" ) {

				print OUT2 "$ped[$i][0] $ped[$i][1] 0 0 $ped[$i][4] $ped[$i][5] ";	#so that "founder's" allele frequencies can be used.

			} else {

				for($k=0;$k<=5;$k++){
					print OUT2 "$ped[$i][$k] ";
				}
			}

			for($c=0;$c<$numchrom;$c++){

				if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

					if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
						next;
					}
				}

				for($k=0;$k<$chr_snp_cnt[$c];$k++){
					#print "Current chromo is ",$c+1,". Snp count is $chr_snp_cnt[$c]\n";
					#print "SNP is $chr_order_snp[$c][$k][0], ind count in geno array is $pedidcnt\n";

					if($whichsamples[$i]!=0){

						$geno=$genos{$chr_order_snp[$c][$k][0]}[$cnt];
						
						if($geno== -1){ 
							print OUT2 "0 0 ";
						}
						if($geno==0){ 
							print OUT2 "1 1 ";
						}
						if($geno==1){ 
							print OUT2 "1 2 ";
						}
						if($geno==2){
							print OUT2 "2 2 ";
						}
					} else {

						print OUT2 "0 0 ";
					}
				}
			}

			if($whichsamples[$i]!=0){	#only increment if previous individual was genotyped!			
				$cnt+=1;
			}
			
			print OUT2 "\n";
		}
	}
	if(defined($plSim)) {
	
		for (my $i = 1; $i <= $plSim; $i++) {

		print "Simulating $i of $plSim individuals for PLINK output.\n";

			print OUT2 "Simulate ";
			print OUT2 "$i ";
			print OUT2 "0 0 "; # pretend that they are all founders...
			print OUT2 "2 0 "; # all females and unaffected

			for($c=0;$c<$numchrom;$c++){

				if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

					if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
						next;
					}
				}

				for($k=0;$k<$chr_snp_cnt[$c];$k++){

#					$snpfreq=$chr_order_snp[$c][$k][4];

					if (rand() <= $chr_order_snp[$c][$k][4]) {
					
						print OUT2 "1 ";

					} else {

						print OUT2 "2 ";
					}

					if (rand() <= $chr_order_snp[$c][$k][4]) {
					
						print OUT2 "1 ";

					} else {

						print OUT2 "2 ";
					}
				}
			}
				
				print OUT2 "\n";
		}
	}

	close(OUT2);
	close(OUT);
}

# output files for personal format, only prints out .pre style files but a) transposed from PLINK and b) per chromosome
# affecteds are clumped together in the final columns, separated by a tab instead of a whitespace from the unaffecteds.
sub print_CP{
my $c;
my $chrom;
my $i;
my $j;
my $k;
my $geno;
my $temp;
my $cnt=0;
my $filename;
my %homoz=();
my $ucnt=0;
	for($c=0;$c<$numchrom;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$filename=$cpDir."/cp_".$chrom.".ped";	
		open(OUT2,">$filename") or print_usage "Cannot write file $filename, $!";
		print OUT2 "Chr\tSNP\trsname\tgenmap\tphysmap\t";
		for($i=0;$i<$pedno;$i++){
			if($whichsamples[$i]!=0 && $ped[$i][5]!=2){ # print out unaffecteds first, then affecteds
				print OUT2 "$ped[$i][1]\t";
			}
		}
		for($i=0;$i<$pedno;$i++){
			if($whichsamples[$i]!=0 && $ped[$i][5]==2){ # print out unaffecteds first, then affecteds
				print OUT2 "$ped[$i][1]\t";
			}
		}
		print OUT2 "S_robdom(D)\tS_HBS(RI)\tS_pairs(R)";
		print OUT2 "\n";
		for($k=0;$k<$chr_snp_cnt[$c];$k++){
			print OUT2 "$chrom\t$chr_order_snp[$c][$k][0]\t"; # chr and snp name (affy style)
			print OUT2 "$chr_order_snp[$c][$k][6]\t"; #rs name
			printf OUT2 "%1.2f", $chr_order_snp[$c][$k][3]; # chr cM pos
			print OUT2 "\t";
			print OUT2 $chr_order_snp[$c][$k][7]; # chr bp pos
			print OUT2 "\t";
			$cnt=0;
			for($i=0;$i<$pedno;$i++){
				if($whichsamples[$i]!=0){
					if($ped[$i][5]!=2){ # print out unaffecteds first, then affecteds
						print OUT2 $genos{$chr_order_snp[$c][$k][0]}[$cnt]+1,"\t";
					}
					$cnt+=1;
				}
			}
			$cnt=0;
			for($i=0;$i<$pedno;$i++){
				if($whichsamples[$i]!=0){
					if($ped[$i][5]==2){ # print out unaffecteds first, then affecteds
						print OUT2 $genos{$chr_order_snp[$c][$k][0]}[$cnt]+1,"\t";
					}
					$cnt+=1;
				}
			}
			# IBS sharing statistics - as adapted from McPeek's proposed IBD sharing statistics		
			# robdom - for robust dominant diseases
			# HBD - for inbred pedigrees
			# S-pairs (second best choice here for recessive diseases in non inbred pedigrees
			# I have tried to rescale all of these by dividing by a maximum so the scores should range between [0,1]
			%homoz=();
			$cnt=0;
			$ucnt=0;
			$homoz{-1}=0;
			$homoz{0}=0;
			$homoz{1}=0;
			$homoz{2}=0;
			for($i=0;$i<$pedno;$i++){
				if($whichsamples[$i]!=0){
					if($ped[$i][5]==2){ 
						$homoz{$genos{$chr_order_snp[$c][$k][0]}[$cnt]}+=1;
						$ucnt+=1;
					}	
					$cnt+=1;
				}
			}
			#robdom
			if(($ucnt-$homoz{-1})==0){
				$temp=-1;
			}
			else{
				#print "7**($homoz{0}+$homoz{1})-1+7**($homoz{2}+$homoz{1})-1 ",7**($homoz{0}+$homoz{1})-1+7**($homoz{2}+$homoz{1})-1,"\n";
				$temp=(7**($homoz{0}+$homoz{1})-1+7**($homoz{2}+$homoz{1})-1)/(2*(7**($ucnt-$homoz{-1}))-2);
			}
			printf OUT2 "%1.4lf",$temp;
			print OUT2 "\t";
			#HBD
			if($homoz{-1}==$ucnt){
				$temp=-1;
			}
			else{
				if($homoz{0}>$homoz{2}){
					$temp=$homoz{0}/($homoz{0}+$homoz{1}+$homoz{2});					
				}
				else{
					$temp=$homoz{2}/($homoz{0}+$homoz{1}+$homoz{2});
				}
			}
			printf OUT2 "%1.4lf",$temp;
			print OUT2 "\t";
			# S_pairs for recessive non inbred pedigrees
			if($homoz{-1}==$ucnt){
				$temp=-1;
			}
			else{
				$temp=(4*($homoz{0}+$homoz{2})+2*(2*$homoz{1}*$homoz{2}+2*$homoz{1}*$homoz{0}+$homoz{1}))/(4*($ucnt-$homoz{-1})**2);
			}
			printf OUT2 "%1.4lf",$temp;
			print OUT2 "\t";
			# all
			# pairs
			print OUT2 "\n";
		}
		close(OUT2);
	}
}

#Output files for Beagle (SR and BL Browning)
sub print_BEAGLE{
my $c;
my $chrom;
my $mfilename;
my $pfilename;
my $j;
my $i;
my $geno;
my $cnt=0;		
	for($c=0;$c<22;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		$pfilename=$beagleDir."/beagle_".$chrom."_".$pfile.".ped";
		$mfilename=$beagleDir."/beagle_".$chrom."_".$pfile.".map";
		open(OUT1, ">$pfilename") or print_usage "Cannot write file $pfilename, $!";
		open(OUT2, ">$mfilename") or print_usage "Cannot write file $mfilename, $!";
		for($j=0;$j<$chr_snp_cnt[$c];$j++){
			print OUT2 "$chr_order_snp[$c][$j][6]\t$chr_order_snp[$c][$j][3]\t1\t2\n"; #rs name/cM/1st allele name/2nd allele name
		}  
		
		close(OUT2);
		print OUT1 "I\tid";
		for($i=0;$i<$pedno;$i++){
			if($whichsamples[$i]!=0){
				print OUT1 "\t".$ped[$i][1]."\t".$ped[$i][1];
			}
		}
		print OUT1 "\n";
		print OUT1 "A\tAffection";
		for($i=0;$i<$pedno;$i++){
			if($whichsamples[$i]!=0){
				print OUT1 "\t".$ped[$i][5]."\t".$ped[$i][5];
			}
		}
		print OUT1 "\n";
		for($j=0;$j<$chr_snp_cnt[$c];$j++){
			$cnt = 0;
			print OUT1 "M\t".$chr_order_snp[$c][$j][0];
			for($i = 0; $i < $pedno; $i++) {
				if($whichsamples[$i]!=0){	
					$geno=$genos{$chr_order_snp[$c][$j][0]}[$cnt];
					if($geno== -1){ 
						print OUT1 "\t0\t0";
					}
					if($geno== 0){ 
						print OUT1 "\t1\t1";
					}
					if($geno== 1){ 
						print OUT1 "\t1\t2";
					}
					if($geno== 2){
						print OUT1 "\t2\t2";
					}	
					$cnt++;
				}
			}
			print OUT1 "\n";
		}
		close(OUT1);
	}		
}


# output files for FEstim 
sub print_FESTIM{
my $c;
my $chrom;
my $i;
my $j;
my $k;
my $geno;
my $temp;
my $snpfreq;
my $cnt=0;
my $filename;
	$filename=$festimDir."/festim_".$pfile."_map";
	open(OUT2,">$filename") or print_usage "Cannot write file $filename, $!";
	# print out the map file for FEstim: This consists of
	# SNP_Name	Chr	Pos(cM) 2	1	Allele1_freq	2	Allele2_freq
	for($c=0;$c<22;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		$chrom=$c+1;
		for($i=0;$i<$chr_snp_cnt[$c];$i++){
			$snpfreq=$chr_order_snp[$c][$i][4];
			if($snpfreq == 0) {
				$snpfreq = 0.001;
			}
			if($snpfreq == 1) {
				$snpfreq = 0.999;
			}
			print OUT2 "$chr_order_snp[$c][$i][6]\t";	#SNP name
			print OUT2 "$chrom\t"; #Chr
			printf OUT2 "%1.4f", $chr_order_snp[$c][$i][3]; #Pos(cM)
			print OUT2 "\t2\t1\t";
			printf OUT2 "%1.3f", $snpfreq;	#Allele1_freq
			print OUT2 "\t2\t";
			printf OUT2 "%1.3f", (1-$snpfreq);	#Allele2_freq
			print OUT2 "\n";
		}
	}
	close(OUT2);
	$filename=$festimDir."/festim_".$pfile."_data";;	
	open(OUT2,">$filename") or print_usage "Cannot write file $filename, $!";
	for($i=0;$i<$pedno;$i++){
		if($whichsamples[$i]!=0){
			for($k=0;$k<=5;$k++){
				print OUT2 "$ped[$i][$k] ";
			}
			for($c=0;$c<22;$c++){

				if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

					if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
						next;
					}
				}

				for($k=0;$k<$chr_snp_cnt[$c];$k++){
					#print "Current chromo is ",$c+1,". Snp count is $chr_snp_cnt[$c]\n";
					#print "SNP is $chr_order_snp[$c][$k][0], ind count in geno array is $pedidcnt\n";
					$geno=$genos{$chr_order_snp[$c][$k][0]}[$cnt];
					if($geno== -1){ 
						print OUT2 "0 0 ";
					}
					if($geno==0){ 
						print OUT2 "1 1 ";
					}
					if($geno==1){ 
						print OUT2 "1 2 ";
					}
					if($geno==2){
						print OUT2 "2 2 ";
					}
				}
			}
			$cnt+=1;
			print OUT2 "\n";
		}
	}
	close(OUT2);
	close(OUT);
}

# printing out RELATE output files. RELATE also takes plink files as input but then uses phys pos instead of genetic map pos for calculations. Hence here we are using RELATE's internal format.
sub print_RELATE{
my $i;
my $j;
my $geno;
my $c;
my $k;
my $cnt=0;
my $filenamerg;
my $filenamerp;
my $filenamerc;
my $filenamero;
	$filenamerg=$relateDir."/relate_".$pfile.".geno";
	$filenamerp=$relateDir."/relate_".$pfile.".pos";
	$filenamerc=$relateDir."/relate_".$pfile.".chr";
	$filenamero=$relateDir."/relate_".$pfile.".opt";
	open(OUTRG,">$filenamerg") || die print "Can't open RELATE genotype file $filenamerg\n";
	open(OUTRC,">$filenamerc") || die print "Can't open RELATE chromosome file $filenamerc\n";
	open(OUTRP,">$filenamerp") || die print "Can't open RELATE positions file $filenamerp\n";
	open(OUTRO,">$filenamero") || die print "Can't open RELATE options file $filenamero\n";
	# print genotype file - data is a matrix of N by M where N=no of inds, M=no of SNPs
	for($i=0;$i<$pedno;$i++){
		if($whichsamples[$i]!=0){
			for($c=0;$c<22;$c++){

			if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c])){    #test if regions are being used and if so whether region included

					if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
						next;
					}
				}

				for($k=0;$k<$chr_snp_cnt[$c];$k++){
					$geno=$genos{$chr_order_snp[$c][$k][0]}[$cnt]+1; # relate wants 0 for missing 1=AA, 2=AB, 3=BB
					print OUTRG "$geno ";
				}
			}
			print OUTRG "\n";
			$cnt+=1;
		}
	}
	# print positions file with all genetic map positions & chromosome file with all 
	for($c=0;$c<22;$c++){

		if ( (%regions && !defined($regions{$c+1})) || !defined($chr_snp_cnt_orig[$c]) ){    #test if regions are being used and if so whether region included

			if ($c <= 21 || ($c == 22 && !defined($regions{"X"}))) {
				next;
			}
		}

		for($k=0;$k<$chr_snp_cnt[$c];$k++){
			print OUTRC $c+1,"\n";
			printf OUTRP "%1.4f", $chr_order_snp[$c][$k][3]; #Pos(cM)
			print OUTRP "\n";
		}
	}
	# option file print out	
	print OUTRO "1 #1=allpairs 0=normal run\n";
	print OUTRO "8 #pair[0]\n";
	print OUTRO "pair[1]\n";
	print OUTRO "0 #double recombination\n";
	print OUTRO "#LD=0=rsq2 LD=1=D //allright everythings gone be allright\n";
	print OUTRO ".025 # min\n";
	print OUTRO "0.001 # alim[0]\n";
	print OUTRO "5.0 # alim[1]\n";
	print OUTRO "0 #doParameter calculation (pars)\n";
	print OUTRO "0.3 # par[0] = a this is only used if doParameter is set to 1\n";
	print OUTRO "0.25 # par[1] = k2 this is only used if doParameter is set to 1\n";
	print OUTRO "0.5 # par[2] = k1 this is only used if doParameter is set to 1\n";
	print OUTRO "16\n";
	print OUTRO "1 #ld_adj\n";
	print OUTRO "0.01 #epsilon\n";
	print OUTRO "50 #back\n";
	print OUTRO " #doPrune\n";
	print OUTRO "0.1 #prune_value\n";
	print OUTRO "0 #fixA\n";
	print OUTRO "0.0 #fixA_value\n";
	print OUTRO "1 #fixK2\n";
	print OUTRO "0.0 #fixk2_value\n";
	print OUTRO "1 #calculateA\n";
	print OUTRO "0.013 #phi_value\n";
	print OUTRO "0.1 #convergence_tolerance\n";
	print OUTRO "3 #times_to_converge\n";
	print OUTRO "8 #times_to_run\n";
	print OUTRO "100 #back2\n";
	print OUTRO "17\n";	
	close(OUTRG);
	close(OUTRC);
	close(OUTRP);
	close(OUTRO);
}

# printing out fastPHASE input files. Scheet and Stephens et al. AJHG 2006
sub print_fastPHASE{
my $i;
my $j;
my $c;
my $p;
my $cnt=0;
my $filenamerg;
my $filenameFPS;
my @ALL;
	$filenameFPS=$fastPHASEDir."/fastPHASE.sh";
	open(OUTFPS,">$filenameFPS") || die print "Can't open fastPHASE shell output file $filenameFPS\n";
	for($c=0;$c<22;$c++){
		$p=$c+1;
		$cnt=0;
		$filenamerg=$fastPHASEDir."/fastPHASE_".$p.".inp";
		open(OUTRG,">$filenamerg") || die print "Can't open fastPHASE input file $filenamerg\n";
		print OUTRG $#which+1;
		print OUTRG "\n";
		print OUTRG "$chr_snp_cnt[$c]\n";
		print OUTRG "P ";
		for($j=0;$j<$chr_snp_cnt[$c];$j++){
			print OUTRG "$chr_order_snp[$c][$j][7] ";
		}
		print OUTRG "\n";
		for($j=0;$j<$chr_snp_cnt[$c];$j++){
			print OUTRG "S";
		}
		print OUTRG "\n";
		for($i=0;$i<$pedno;$i++){
			if($whichsamples[$i]!=0){
				print OUTRG "# $ped[$i][1]\n";
				for($j=0;$j<$chr_snp_cnt[$c];$j++){
					if($genos{$chr_order_snp[$c][$j][0]}[$cnt]== -1){
						$ALL[$j][0]="?";
						$ALL[$j][1]="?";
						next;
					}
					if($genos{$chr_order_snp[$c][$j][0]}[$cnt]==0){
						$ALL[$j][0]=1;
						$ALL[$j][1]=1;
						next;
					}
					if($genos{$chr_order_snp[$c][$j][0]}[$cnt]==1){
						$ALL[$j][0]=1;
						$ALL[$j][1]=2;
						next;
					}
					if($genos{$chr_order_snp[$c][$j][0]}[$cnt]==2){
						$ALL[$j][0]=2;
						$ALL[$j][1]=2;
						next;
					}
				}
				for($j=0;$j<$chr_snp_cnt[$c];$j++){
					print OUTRG "$ALL[$j][0]";
				}
				print OUTRG "\n";
				for($j=0;$j<$chr_snp_cnt[$c];$j++){
					print OUTRG "$ALL[$j][1]";
				}
				print OUTRG "\n";
				$cnt+=1;
			}						
		}
		@ALL=();
		# write vanilla shell scrip to run fastPHASE.
		print OUTFPS "fastPHASE -ochr".$p." fastPHASE_".$p.".inp\n";
	}
}
			

# prints out the settings that were used for this LINKDATAGEN run
sub config {
	print "\nConfiguration of LINKDATAGEN:\n";
	
	print "\nMandatory options:\n\n";
	
	if (defined($annotFile)) {
		print "-annotFile\t$annotFile\n";
	}
	if (defined($chip)) {
		print "-chip\t\t\t$chip\n";
	}
	if (defined($callDir)) {
		print "-callDir\t\t$callDir\n";
	}
	if (defined($callFile)) {
		print "-callFile\t\t$callFile\n";
	}
	if (defined($data)) {
		print "-data\t\t\t$data\n";
	}
	if (defined($pedfile)) {
		print "-pedfile\t\t$pedfile\n";
	}
	if (@prog) {
		print "-prog\t\t\t@prog\n";
	}
	if (defined($popHetTest) && !defined($bestPopTest)) {
		print "-popHetTest\t\t$popHetTest\n";
	}
	if (defined($freq)) {
		print "-freq\t\t\tSELECTED\n";
	}
	if (defined($whichSamplesFile)) {
		print "-whichSamplesFile\t$whichSamplesFile\n";
	}
	if (defined($whichSamplesList)) {
		print "-whichSamplesList\t$whichSamplesList\n";
	}
	if (defined($MEO)) {
		print "-MendelErrors\t\t$MEO - ";

		if ($MEO eq "keep") {
			print "Mendelian errors are only reported.  They have not been removed from your dataset.\n";

		} elsif ($MEO eq "removeSNPs") {
			print "any SNP displaying a Mendelian error has been removed completely from your dataset - across all families.\n";

		} elsif ($MEO eq "removeGPF") {
			print "genotypes for any SNP displaying a Mendelian error have been wiped (\"made missing\") on a per-family basis (i.e. for only the family displaying the Mendelian error).\n"; 
		}
	}
	
	print "\nOthers:\n\n";
	
	if (defined($actg)) {
		print "-actg\t\t\t$actg\n";
	}
	if (defined($annotDir)) {
		print "-annotDir\t\t$annotDir\n";
	}
	if (defined($bestPopTest)) {
		print "-bestPopTest\t\tTest for best population performed.  See results in table above.\n";
	}
	if (defined($binsize)) {
		print "-binsize\t\t$binsize cM\n";
	}
	if (defined($crlmm)) {
		print "-crlmm\t\t\t$crlmm\n";
	}
	if (defined($fileKeepSNPs)) {
		print "-fileKeepSNPs\t\t$fileKeepSNPs\n";
	}
	if (defined($fileRemoveSNPs)) {
		print "-fileRemoveSNPs\t\t$fileRemoveSNPs\n";
	}
	if (defined($help)) {
		print "-help\t\t\t$help\n";
	}
#	if (defined($keepME)) {
#		print "-keepME\t\t\tSELECTED - SNP markers with simple Mendelian errors have been kept (not removed) from your dataset.\n";
#	}
#	else {
#		print "-keepME\t\t\tNOT SELECTED - instead, SNP markers displaying simple Mendelian errors have been removed.\n";
#	}
	if (defined($minDist)) {
		print "-minDist\t\t$minDist cM";
		if($minDist <= 0.15) {
			print " - WARNING: Minimum distance is very small. Linkage disequilibrium may be a problem for you.";
		}
		elsif($minDist > $binsize) {
			print " - WARNING: Minimum distance > binsize.  This is counter intuitive, but maybe what you want.";
		}
		print "\n";
	}
	if (defined($noX)) {
		print "-noX\t\t\t$noX\n";
	}
	if (defined($outputDir)) {
		print "-outputDir\t\t$outputDir\n";
	}
	
	if (!defined ($bestPopTest)) {
	
		if (defined($pop)) {
			print "-pop\t\t\t$pop\n";
		}
		if (defined($popCol)) {
			print "-popCol\t\t\t", ($popCol + 1), "\n";
		}
	}
	
	if (defined($randomSNP)) {
		print "-randomSNP\t\tSELECTED - a random SNP will be selected from each bin interval.\n";
	}
	else {
		print "-randomSNP\t\tNOT SELECTED - instead, the most heterozygous SNP will be selected from each bin interval.\n";
	}
	if (defined($regions)) {
		print "-regions\t\t$regions\n";
	}
	if (defined($regionsFile)) {
		print "-regionsFile\t\t$regionsFile\n";
	}
	if (defined($removeAIS)) {
		print "-removeAIS\t\tSELECTED - SNP markers absent from any Illumina input file have been removed.\n";
	}
	if (defined($removeWFHBS)) {
		print "-removeWFHBS\t\t$removeWFHBS";
	
		if ($removeWFHBS eq "u") {
			print " - union (across families) of within-family homozgyosity-by-state SNP markers will be removed.\n";
		}
		elsif ($removeWFHBS eq "i") {
			print " - intersection (across families) of within-family homozgyosity-by-state SNP markers will be removed.\n";
		}
	}
	else {
		print "-removeWFHBS\t\tNOT SELECTED - within-family homozgyosity-by-state SNP markers will NOT be removed.\n";
	}
	if (defined($seed)) {
		print "-seed\t\t\t$seed\n";
	}
	if (defined($minMAF)) {
		print "-minMAF\t\t\t$minMAF (note - within LINKDATAGEN this is also equivalent to " . (1 - $minMAF) . " which is 1 - $minMAF).\n";
	}
	if (defined($maxMAF)) {
		print "-maxMAF\t\t\t$maxMAF (note - within LINKDATAGEN this is also equivalent to " . (1 - $maxMAF) . " which is 1 - $maxMAF).\n";
	}
	if (defined($plSim)) {
		print "-plSim\t\t\t$plSim individual(s) will be simulated in the PLINK dataset given the selected population allele frequencies.\n";
	}
	if ($data eq "i") {
	
		if(defined($kCT)) {
			print "-kCT\t\t\tSELECTED - complementary transverion SNP markers (A/T, T/A, C/G, G/C) were not excluded from the Illumina annotation file.\n";
		}
		else {
			print "-kCT\t\t\tNOT SELECTED - complementary transverion SNP markers (A/T, T/A, C/G, G/C) have been excluded from the Illumina annotation file.\n";
		}
	}
}

