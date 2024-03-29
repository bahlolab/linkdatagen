#!/usr/bin/perl -w
# svn $Revision: 997 $
# svn $LastChangedDate: 2016-02-03 13:38:44 +1100 (Wed, 03 Feb 2016) $

#This is a companion script to linkdatagen.pl, to be used on .vcf or .vcf-like files.
#This script creates a brlmm genotype call file from genotypes called from MPS data
#BRLMM call format:
#-1 = No call
#0 = AA call
#1 = AB call
#2 = BB call

=head1 VCF2LINKDATAGEN


                    ----------VCF2LINKDATAGEN----------
                           CJ Bromhead, M Bahlo
                    -----------------------------------

=head1 Usage

  vcf2linkdatagen.pl -idlist <idlist> -annotfile <annotfile> [-pop -missingness -mindepth -min_MQ -min_FQ -minP_strandbias -minP_baseQbias -minP_mapQbias -minP_enddistbias -alleleACGT]

=head2 Required parameters

  -idlist <idlist> - <idlist> is a file containing a list of filenames & paths to input VCF files

  -annotfile <annotfile> - <annotfile> is the file containing the annotation data

  -variantCaller / -vc <mpileup | unifiedGenotyper | mp | ug> - The variant caller used in genotype calling, where
    SAMtools mpileup use mpileup or mp (at present this only supports the older SAMtools 0.1.19, not 1.0.0+)
    GATK UnifiedGenotyper use unifiedGenotyper or ug

=head2 Optional parameters: [default value in square brackets]

  -pop [CEU] - Population frequency choice. Default is CEU. Choose one of:
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
    ALL : No filtering performed on population allele frequencies.

  -missingness [1] the maximum proportion of missing genotype calls for a SNP to be output to the brlmm file.

  -extra_info : If defined, vcf2linkdatagen will print out a file for each VCF, vcf_ID_callswithextrainfo.  This file contains the fields rs_name, genocall, chr, pos, DP, sumDP4, MQ, FQ and AF1 for each called genotype

=head2 Thresholds for genotype calling: [default value in square brackets]

A vcf genotype call can be converted into a brlmm genotype call provided that the values in the information field of the vcf file do not fall below the following thresholds:

    -mindepth [10] An integer value for the minimum read depth at a site.  This is calculated from the DP4 field in the vcf file rather than the DP field, as the DP4 field counts only high quality base calls.
    -min_MQ [10] The minimum root mean squared mapping quality (MQ) at a site.
    -min_FQ [10] The minimum absolute value of consensus quality (FQ) at a site.

The following four options are thresholds for values of PV4.  Not all lines of the vcf file will contain a PV4 field but for those that do we can set thresholds.  The PV4 values are P-values for 1) strand bias (exact test); 2) baseQ bias (t-test); 3) mapQ bias (t-test); 4) tail distance bias (t-test).

    -minP_strandbias [0.0001] minimum P value for strand bias at a site.
    -minP_baseQbias [1e-100] minimum P value for base Q bias at a site.
    -minP_mapQbias [0] minimum P value for map Q bias at a site.
    -minP_enddistbias [0.0001] minimum P value for tail distance bias at a site.

=cut 

use 5.006;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $script_revision = "1000";
'$LastChangedDate: 2020-08-06 (Sat, 06 Aug 2022) $' =~ /\((.+)\)/;
my $last_changed_date = $1;

# Compressed file opening modules, allows script to work even if modules are not installed:
my $perliogzip = 0;
eval { require PerlIO::gzip };
unless($@) { $perliogzip = 1 };
my $perlioviabzip2 = 0;
eval { require PerlIO::via::Bzip2 };
unless($@) { $perlioviabzip2 = 1 };

my $annot_dir;
my $annotfile;
my $min_MQ;
my $mindepth;
my $DP;
my $min_FQ;
my $minP_strandbias;
my $minP_baseQbias;
my $minP_mapQbias;
my $minP_enddistbias;
my $min_GQ;
my $idlist;
my $log;
my $pop;
my $missingness;  	#new as of May 5, this means maximum proportion of missing genotype calls allowed for a 			
						#SNP to be output to the brlmm.  Must be between 0 and 1.  Default to 1 implies that at 			
						#the filtering stage vcf2linkdatagen will not discard any calls.  A missingness value of 				
						#0.5 would mean that any SNP with >50% missing genotype data would not be included in the 				
						#output. *LOG
my $allele_ACGT; #not working yet, July 2011
my $help;
my $extra_info;
my $variantCaller;
my $vc;
my $copen_gz_warn;

GetOptions(
	"min_MQ=i"=>\$min_MQ,
	"mindepth=i"=>\$mindepth,
	"annotfile=s"=>\$annotfile,
	"minP_strandbias=f"=>\$minP_strandbias,
	"minP_baseQbias=f"=>\$minP_baseQbias,
	"minP_mapQbias=f"=>\$minP_mapQbias,
	"minP_enddistbias=f"=>\$minP_enddistbias,
	"min_FQ=f"=>\$min_FQ,
	"min_GQ=f"=>\$min_GQ,
	"log"=>\$log,
	"missingness=f"=>\$missingness,
	"annot_dir=s"=>\$annot_dir,
	"idlist=s"=>\$idlist,
	"allele_ACGT"=>\$allele_ACGT,
	"help"=>\$help,
	"pop=s"=>\$pop,
	"extra_info"=>\$extra_info,
	"variantCaller=s"=>\$variantCaller, 
	"vc=s"=>\$vc ) or
	print_usage("Error in parsing options: unknown option/s given (see top of this message for more details).");

#fill in undefined values
if(!defined($annotfile) ) 		{$annotfile="annotHapMap2.txt";}
if(defined($annot_dir) ) 		{$annotfile="$annot_dir/$annotfile";}# add directory to annotfile if annot_dir given
if(!defined($min_MQ) ) 		{$min_MQ=10;}
if(!defined($mindepth) ) 	{$mindepth=10;}
if(!defined($min_FQ) )		{$min_FQ=10;}
if(!defined($min_GQ) )		{$min_GQ=20;}
if(!defined($minP_strandbias) )	{$minP_strandbias=0.0001;}
if(!defined($minP_baseQbias) )	{$minP_baseQbias=1e-100;}
if(!defined($minP_mapQbias) )	{$minP_mapQbias=0;}
if(!defined($minP_enddistbias))	{$minP_enddistbias=0.0001;}
if(!defined($missingness) )		{$missingness=1;}

if(defined($help)){
    pod2usage(
        -verbose => 2, 
        -width => 50,
        -noperldoc => ! -t STDOUT,
    );
}

print STDERR "\n\t\t----------VCF2LINKDATAGEN----------\n\t\t\tCJ Bromhead, M Bahlo\n\t\t-----------------------------------\n\n";
print STDERR "This script updated on $last_changed_date (r$script_revision)\n";

if(!defined($annotfile)) {
	print_usage("Missing annotation file");
}

if(!defined($idlist) && !defined($ARGV[0])) {
	print_usage("Missing idlist or annotation file");
}

if(@ARGV > 1) {
    print_usage("Only one VCF file may be given directly. Use the idlist option instead.\nThis error may be due to missing hyphens in other options or options without parameters");
}

if(!(defined($variantCaller) || defined($vc))) {

	print_usage("You must specify either -variantCaller [mpileup | unifiedGenotyper] or in short form -vc [mp | ug]");
}

if (defined($variantCaller) && defined($vc) && $variantCaller != $vc) {

	print_usage("You have speficied both -variantCaller and -vc but with different values.  You only need specify -variantCaller or -vc.");
}

if (defined($vc) && !defined($variantCaller)) {

	$variantCaller = $vc;
}

$variantCaller = lc $variantCaller;

if (!($variantCaller eq "unifiedgenotyper" || $variantCaller eq "ug" || $variantCaller eq "mpileup" || $variantCaller eq "mp")) {

	print_usage("Must specify either -variantCaller [mpileup | unifiedGenotyper] or in short form -vc [mp | ug]");

} elsif ($variantCaller eq "ug") {

	$variantCaller = "unifiedgenotyper";
	
} elsif ($variantCaller eq "mp") {

	$variantCaller = "mpileup";
}

my %flipper=();
$flipper{"A"}="T";
$flipper{"T"}="A";
$flipper{"G"}="C";
$flipper{"C"}="G";
$flipper{"N"}="N";

my %genos_orig=();
my %annot_orig=();
my $num_files;
sub filter_missing_calls_vcf;

if(defined($log)) {	
	my $logfile = "vcf2linkdatagen.log";
	copen(*LOG, ">$logfile") or die "Can't open $logfile for writing. $!";
}

open(*OUT, ">-") or die "$!";

# -------------- populations -------------------
my %pops=("CEU",1,"ASW",2,"CHB",3,"CHD",4,"GIH",5,"JPT",6,"LWK",7,"MEX",8,"MKK",9,"TSI",10,"YRI",11);
my @hapmap2=(1,3,6,11);

my $tmp1="annotHapMap2";
my $tmp2="annotHapMap3";
my $popcol; 
my $key;

print STDERR "Annotation file $annotfile will be used to specify which SNPs go into the brlmm.out file\n\n";

if( !defined($pop) ) {
	print STDERR "No population has been specified for the frequencies: Caucasians (CEU) chosen by default\n";
	$pop="CEU";
	$popcol="CEU"; 
} elsif ($pop eq "ALL") {

	$popcol="CEU";	#when pop=ALL, it doesn't matter what popcol number is really, just anything to keep the code working

	print STDERR "\"ALL\" populations have been selected.  This means SNPs will not be filtered on the basis of allele frequencies for any population.\n";

} else {
	print STDERR "Population $pop has been specified for the frequencies.\n";

	if(defined($annotfile) && $annotfile=~/$tmp1/){
		foreach $key (@hapmap2){

			if($pops{$pop}==$key){ # is one of the hapmap2 pops and ok
				goto POPCOL;
			}
		}

		print_usage("For MPS data currently the only choice is the CEU, CHB, JPT or YRI with HapMap2 allele frequencies.\n");
	}

	if(defined($annotfile) && $annotfile=~/$tmp2/){

		foreach $key (@hapmap2){

			if($pops{$pop}==$key){ # is one of the hapmap2 pops
				print STDERR "WARNING: Use annotHapMap2.txt annotation file instead, for more SNPs.\n";
			}
		}
	}

	POPCOL: $popcol = $pop; 
}

copen(*ANNOT, '<', $annotfile) || die "Can't open annotation file $annotfile. $!";

#VCF2LINKDATAGEN.PL MAIN
#------------------------------------------------------------------------------------
print STDERR "Started at ";
&print_time();
print STDERR "\n\n";
#First we must read in the annotation file
&read_in_annot_vcf();
print STDERR "-----------------------------------------------------------------\n";

#Then we read in the genotype data
&read_in_vcf();
print STDERR "-----------------------------------------------------------------\n";

#then a recoding of the annot hash and genotype hash to make the key value the rs ID of the SNP
#in conjunction with discarding SNPs from the annot hash that are not in the geno hash
&recode_hashes_vcf();
print STDERR "-----------------------------------------------------------------\n";
print STDERR "Recoding genotypes to BRLMM format\n";
&recode_geno_vcf();	
print STDERR "-----------------------------------------------------------------\n";

#Now we must filter calls based on missingness threshold if $missingness < 1
if($missingness < 1) {
	#print STDERR "**7\n";
	#print "Missingness filter:\n";
	filter_missing_calls_vcf();
	print STDERR "-----------------------------------------------------------------\n";
}

#and finally we write the brlmm file
&write_brlmm();		
print STDERR "\nFinished at ";
&print_time();
print STDERR "\n";	

#-------------------------------------------------------------------------------------

#sub routines
sub print_usage {
    die("\nVCF2LINKDATAGEN.pl has aborted.\n\n$_[0]\n\nUse -help for further instructions, or refer to the LINKDATAGEN manual.\n");
}

sub array_index{   #new subroutine (May 2011)
	my $c;
	my %index_of_array1=();
	for($c = 0; $c < @_; $c++) {	#here @_ is the array we wish to index, and $_ is used to refer directly
		$index_of_array1{$_[$c]} = $c;			#to elements of the array
	}
	return %index_of_array1;
}

# read in annotation files are now taken from our internal format as generated by Katherine Smith, Rick Tankard and Catherine Bromhead
#Chrom	
#physical_position_(bp)	
#rs_name	
#Strand	
#TOPBOT strand
#deCODE_genetic_map_position	
#alleleA/alleleB
#allele_frequencies_1)CEU 2)ASW	3)CHB 4)CHD 5) GIH 6)JPT 7)LWK 8)MEX 9)MKK10)TSI 10)YRI

# Chrom physical_position_build37 rs_name Strand TOPBOT deCODE_genetic_map_position A B allele_frequencies_CEU ASW CHB CHD GIH JPT LWK MEX MKK TSI YRI

sub read_in_annot_vcf(){
	print STDERR "***********READ IN ANNOTATION SUBROUTINE******************\n";
	my $i;
	my $c;
	my $line;
	my $line_cnt=-1;
	my $key;
	my $test=0;
	my $rsname_col; 
	my $chr_col;
	my $strand_col;
	my $genetic_pos_col;
	my $physical_pos_col;
	my $freq_col;
	my $alleleA_col;
	my $alleleB_col;
	my $earlier_annotfile=1;
	my $key1;
	my @temp;
	my @arr4;
	my $totalline_cnt=0;
	my %hash_header_index=();
	my @header_array;
    my $head_comments = 1;
    my $annot_version_check = 0;
	ANNOTLOOP: while(<ANNOT>){
        chomp;
		$totalline_cnt+=1;
        if($head_comments) {
            # Check head comment lines
            if(/^#/) {
                # In comment header
                if(/^# Minimum LINKDATAGEN revision: (\d+)$/) {
                    if($script_revision < $1) {
                        # this script is too old
                        print_usage("Annotation file is requires newer revision of LINKDATAGEN. (Revision $1 required, but have $script_revision)"); 
                    }
                    $annot_version_check = 1;
                }
                next ANNOTLOOP;
            } else {
                $head_comments = 0;
                unless($annot_version_check) {
                    print_usage("Annotation minimum version not found. Are you using the latest annotation file?");
                }
            }
        }
		@temp=split(/\s+/,$_);
		if($line_cnt == -1) {
			@header_array = @temp;
			%hash_header_index = array_index(@header_array);
            $rsname_col=$hash_header_index{"rs_name"};
			$chr_col=$hash_header_index{"Chrom"};
			$strand_col=$hash_header_index{"Strand"};
			$genetic_pos_col=$hash_header_index{"deCODE_genetic_map_position"};
			$physical_pos_col=$hash_header_index{"physical_position_build37"} || $hash_header_index{"physical_position_build38"};
            $freq_col = defined($hash_header_index{"$popcol"}) ? $hash_header_index{"$popcol"} : $hash_header_index{"allele_frequencies_$popcol"}; 
            $alleleA_col = $hash_header_index{"A"};
            $alleleB_col = $hash_header_index{"B"};
			if(defined($rsname_col) &&  defined($chr_col) &&  defined($strand_col) &&  defined($genetic_pos_col) &&  defined($physical_pos_col) &&  
                defined($freq_col) ){
				print STDERR "Matched column indices in annotation file\n";
			}
			else{
				print_usage("Unable to match column indices in annotation file. The file may be corrupted.");
			}
            $line_cnt+=1;
            next ANNOTLOOP;
		}
		if( ( $temp[$freq_col] =~/\d+/ ) || ( $pop eq "ALL" ) ){	
		#if($temp[$freq_col] =~/\d+/ || $temp[$freq_col] ne "NA"){	# I don't know why this is eq, changed to ne - MBjan2012		
			$key1 = $temp[$chr_col]."pos".$temp[$physical_pos_col];  #key for hash will change
			$annot_orig{$key1}[4]=$temp[$freq_col]; #May 5 pops option will come later
            $annot_orig{$key1}[2] = $temp[$alleleA_col]."/".$temp[$alleleB_col];
			$annot_orig{$key1}[1]=$temp[$chr_col]; # chromosome # TODO: check that now having the "chr" in front doesn't affect this script
			$annot_orig{$key1}[0]=$temp[$rsname_col]; # snp name
			$line_cnt+=1;	
			if($line_cnt%100000==0){
				print STDERR "Read in $line_cnt SNP annotation lines\n";
			}	
		}
	}
	print STDERR "# of SNPs in the annotation file = $totalline_cnt\n";
	print STDERR "# of SNPs with allele frequency data for population $pop in annotation file = $line_cnt\n"; 
}

sub read_in_vcf() {
	
#Now we read through our vcf file
#chr1	888659	.	T	C	226	.	DP=26;AF1=1;CI95=1,1;DP4=0,0,9,17;MQ=49;FQ=-81	GT:PL:GQ	1/1:234,78,0:99
#chr1	990380	.	C	.	44.4	.	DP=13;AF1=7.924e-09;CI95=1.5,0;DP4=3,10,0,0;MQ=49;FQ=-42	PL	0

	#we see whether idlist is defined.  If it is, we read the list of vcf file names from idlist
	my @vcf_ids=();
	my @arr4;
	my $extra_info_file;
	my $c;
	if(defined($idlist)) {
		copen(*ID,"<$idlist") || die "Can't open idlist file $idlist. $!";
		@arr4=<ID>;
		print STDERR "Reading in the idlist file...\n";
		@vcf_ids = map {
			chomp;
			split(/\r/); # allows Mac-style line breaks
		} @arr4;
		print STDERR join("\n", @vcf_ids) . "\n";
		$num_files = @vcf_ids;
		undef(@arr4);
		close(ID);
	}
	else{
		$vcf_ids[0] = "$ARGV[0]";
		$num_files = 1;
	}
	print STDERR "\n";
	for($c = 0; $c < @vcf_ids; $c++) {
		copen(*IN, '<', $vcf_ids[$c]) or die "Can't open $vcf_ids[$c]. $!";
		print STDERR "Reading in VCF file ".$vcf_ids[$c]."\n";
		if(defined($extra_info)) {
			$extra_info_file = $vcf_ids[$c]."_callswithextrainfo\n";
			copen(*EXTRA, ">$extra_info_file") or die "Can't open $extra_info_file for writing. $!";
			print EXTRA "rs_name\tgenocall\tchr\tpos\tDP\tsumDP4\tMQ\tFQ\tAF1\n";
		}
		my $variantCallerDetected = 0; 
		while(<IN>) {
			if( $_ =~ /^#/ ){
				if(/^##samtoolsVersion=/) {
					$variantCaller eq "mpileup" or die "VCF appears to be produced by SAMtools, not $variantCaller as specified.\n";
					$variantCallerDetected = 1;
				}
				if(/^##GATKCommandLine=<ID=UnifiedGenotyper,/) {
					$variantCaller eq "unifiedgenotyper" or die "VCF appears to be produced by UnifiedGenotyper, not $variantCaller as specified.\n";
					$variantCallerDetected = 1;
				}
				if(/^##GATKCommandLine=<ID=HaplotypeCaller,/) {
					$variantCaller eq "haplotypecaller" or die "VCF appears to be produced by HaplotypeCaller, not $variantCaller as specified.\n";
					$variantCallerDetected = 1;
				}
				if ( /^#[^#]/ ) {
					unless($variantCallerDetected) {
						warn "Variant detection algorithm was not found in the VCF header!\n";
					}
				}
				next;
			}
			my $ind_counter=0;
			my $skip = 0;
			my @tmp = split(/\s+/,$_);		
			#step one: is this a SNP in our annotation file?  Otherwise we can skip it.
			my $chr = $1 if ($tmp[0] =~ /chr([\S]+)/);
			my $pos = $tmp[1];
			my $key1 = "chr".$chr."pos".$pos;  #we use this value to locate the SNP if it is in our annotation file
			if(!(defined($annot_orig{$key1}[1]))) {
				#print STDERR "skipping ".$key1."\n";
				$skip = 1;
			}
			#else{ print STDERR "recognised key ".$key1."\n"; }
			if($skip == 0) {  #do not continue if we have decided to skip this line
				my $totaldepth = 0;
				my $MQ = 0;
				my $rearranged = 0;  #my variable to keep track of whether we rearrange the ref and alt alleles
				my $brlmm;
				my $geno;
				#Now we have $tmp[0]=CHROM, $tmp[1]=POS, $tmp[2]=ID, $tmp[3]=REF, $tmp[4]=ALT, $tmp[5]=QUAL, $tmp[6]=FILTER, $tmp[7]=INFO, $tmp[8]=FORMAT
				$pos = $tmp[1];
				my $ref = $tmp[3];
				my $alt = $tmp[4];
				my $qual = $tmp[5];
				my $info = $tmp[7];
				my $format = $tmp[8];
				my $sample_1 = $tmp[9];
				#Now we grab bits from $info using pattern matching!

				if ($variantCaller eq "mpileup") {

					$DP = $1 if ($info =~ /DP=([\d]+)/);
				}

				my $FQ = $1 if ($info =~ /FQ=(-?[\d\.]+)/);
				$MQ = $1 if ($info =~ /MQ=([\d]+)/);
				my $AF1 = $1 if ($info =~ /AF1=([^;]*);/);
				my @PV4 = ($1, $2, $3, $4) if ($info =~ /PV4=([^\,]+),([^\,]+),([^\,]+),(\S+)/);
				my @DP4 = ($1, $2, $3, $4) if ($info =~ /DP4=(\d+),(\d+),(\d+),(\d+)/);
				if(defined($DP4[0])) {
					$totaldepth = $DP4[0] + $DP4[1] + $DP4[2] + $DP4[3];
				}
				# Extracting format field
				my @formats = split(/:/, $format);
				my @formatinfo_1 = split(/:/, $sample_1);
				if(@formats != @formatinfo_1) { die "Sample info has different number of fields from sample data. FORMAT: $format INFO: $sample_1" };
				my %format_hash; 
				@format_hash{@formats} = @formatinfo_1;
							
				#NOW TO CHECK THAT THINGS ARE TURNING OUT OK
				if(defined($log)) {
					print LOG substr($_, 0, length($_) - 1)."\t";
					#print LOG "chr = ".$chr.", FQ = ".$FQ.", MQ = ".$MQ.", AF1 = ".$AF1.", DP4 = ".$DP4[0].",".$DP4[1].",".$DP4[2].",".$DP4[3].",";
					#if(defined($PV4[0])) {
						#print LOG "PV4 = ".$PV4[0].",".$PV4[1].",".$PV4[2].",".$PV4[3];
						#}
					#print LOG "\t";
					#extra check for FQ
					#	print LOG "FQ  ".$FQ." ";
				}
				#if(abs($FQ) < $min_FQ) {
					#print LOG "REJECT\t";
					#}
				#else{
					#print LOG "ACCEPT\t";
					#}
				$key1 = "chr".$chr."pos".$pos;  #since key is already taken - we use this value to locate the SNP if it is in our annotation file

				if ($variantCaller eq "mpileup") {
					if( !defined($FQ)
						|| ($tmp[7] =~ /INDEL/ || $MQ < $min_MQ || $totaldepth < $mindepth || length($alt) > 1 || abs($FQ) < $min_FQ)
						|| (defined($format_hash{GQ}) && $format_hash{GQ} < $min_GQ)) {
							$skip = 1;
					} #5 conditions for skipping: If it is an indel we skip, 
				}
				elsif ($variantCaller eq "unifiedgenotyper") {
				
					if (defined ($format_hash{"GT"}) && defined ($format_hash{"DP"})) {

						if ($format_hash{"GT"} eq "./.") {

							$skip = 1;

						} else {

						$DP = $format_hash{"DP"};

							if( 
								($tmp[7] =~ /INDEL/ || $MQ < $min_MQ || $DP < $mindepth || length($alt) > 1 )
								|| (defined($format_hash{GQ}) && $format_hash{GQ} < $min_GQ)) {
									$skip = 1;
							} #5 conditions for skipping: If it is an indel we skip, 
						}
					} else {

						$skip = 1;					
					}
					
				}

				#if the quality, depth or abs(FQ) is below our threshold we skip and if there is more than one alternate allele we skip.
				if(defined($log)) {
					if($skip == 0) {print LOG "ACCEPT\t"; }
					if($skip == 1) {print LOG "REJECT\t"; }
					if($tmp[7] =~ /INDEL/) {print LOG "INDEL*\t"; }
					if(!(defined($annot_orig{$key1}[1]))) { print LOG "!ANNOT*\t"; }
					#if((defined($annot_orig{$key1}[1]))) { print LOG "!#####***\t".$key1."\t"; }
					if($variantCaller eq "mpileup" && $totaldepth < $mindepth)  {print LOG "DEPTH*\t"; }
					if($variantCaller eq "unifiedgenotyper" && $DP < $mindepth) {print LOG "DEPTH*\t"; }
					if(length($alt) > 1) {print LOG "ALT*\t"; }
					if(!defined($FQ)) {print LOG "FQ=nan\t"; } elsif (abs($FQ) < $min_FQ) {print LOG "FQ*\t"; }
					if(defined($format_hash{GQ}) && $format_hash{GQ} < $min_GQ) {print LOG "GQ*\t"; }
					if($MQ < $min_MQ) {print LOG "MQ*\t"; }
					print LOG "\n";
				}
				if($skip == 0) {
					if(defined($PV4[3])) {	
						if($PV4[0] < $minP_strandbias || $PV4[1] < $minP_baseQbias || $PV4[2] < $minP_mapQbias || $PV4[3] < $minP_enddistbias) {
							$skip = 1;
						}
					}  #If we have PV4 values defined we check that they are not below our thresholds.  If they are below, we skip.		
				}
				if($skip == 0) {
					if(defined($format_hash{GT}) && $format_hash{GT} =~ /\.\/\./) {
						$skip = 1;
					}
				}
				if($skip ==0 ) {	 #only continue if we have not decided to skip this SNV
					my $annot_allA = $annot_orig{$key1}[7+11]; #extract reference alleles from annotation file
					my $annot_allB = $annot_orig{$key1}[8+11];
					if($variantCaller eq "mpileup" && $FQ == 0) {
						$skip = 1;
						#print "Skipping genotype at chr ".$chr." position ".$pos." due to FQ=0\n";
					}
					#now assign a genotype for the called variant, based on the GT field
					# Only SNV sites with one alternative should be present here if previous filters have worked
					if(defined($format_hash{GT})) {
						$format_hash{GT} =~ /^([01])[\/|]([01])$/ && ($1 < 2) && ($2 < 2) or die "Bad GT field $format_hash{GT}.";
						$geno = ( $1 ? $alt : $ref ) . ( $2 ? $alt : $ref ); 
					} else {
						$geno = $ref.$ref;
					}

					if(defined($geno)){				
						if(ord(substr($geno, 0, 1)) > ord(substr($geno, 1, 1))) {
							$geno = substr($geno, 1, 1).substr($geno, 0, 1);  #put the hets in alphabetical order
						}
						$genos_orig{$key1}[$c] = $geno;  #will be recoded to brlmm in another subroutine
						#print STDERR "assigning geno for person ".$c.", key ".$key1.":  ".$geno."\n";
						if(defined($log)){
							print LOG $geno."\n";
						}
						undef($geno);  #for next loopy
					}
					#now would be a good time to print to the "katfile";
					if(defined($extra_info)) {
						print EXTRA $annot_orig{$key1}[0]; # print rs name to Katherine's special file
						if(defined($genos_orig{$key1}[$c])) {
							print EXTRA "\t".$genos_orig{$key1}[$c];
						}
						else{
							print EXTRA "\t-1";
						}
						if($variantCaller eq "mpileup") {print EXTRA "\t".$chr."\t".$pos."\t".$totaldepth."\t".$DP."\t".$MQ."\t".$FQ."\t".$AF1."\n"};
						if($variantCaller eq "unifiedgenotyper") {print EXTRA "\t".$chr."\t".$pos."\t".$DP."\t".$DP."\t".$MQ."\t".$FQ."\t".$AF1."\n"};
					}
				}
			}
		}
	}
}
	
sub write_brlmm(){
	print STDERR "Writing brlmm file\n";	
	my $c;
	my $key;
	my @callarray=();
	foreach $key (keys %genos_orig){
		print OUT $annot_orig{$key}[0];
		for($c = 0; $c < $num_files; $c++) {  #c<1 for now, multiple people to come
			print OUT "\t".$genos_orig{$key}[$c];
		}
		print OUT "\n";
	}
}

sub recode_hashes_vcf() {
my $key;  #old key chr4pos12345678
my $newkey;  #new key rs98777777
my $genocalls;  #default to 0, set to 1 if SNP has genotype calls
my $c;
my %temp_hash_annot=();
my %temp_hash_geno=();
my $snp_miss_cnt=0;
my $key_cnt=0;
my $snps_with_data=0;
	print STDERR "In the recoding hash subroutine\n"; #ERROR CHECKING
	foreach $key (keys %annot_orig){
		$genocalls = 0;
		for($c=0;$c<$num_files;$c++){
			if(defined($genos_orig{$key}[$c])){
				$genocalls = 1;
			}
		}
		if($genocalls==0) {
			#print STDERR "we found that this SNP is in the annot file but not in the geno file\t"; #ERROR CHECKING
			#print STDERR $annot_orig{$key}[0]."\n";
			#delete($annot_orig{$key});
			$snp_miss_cnt+=1;
		}
		if($genocalls == 1){
			#print STDERR "looplooploop\n"; #ERROR CHECKING
			$newkey = $annot_orig{$key}[0];
			#print STDERR "new key! ".$newkey."\n";
			$temp_hash_annot{$newkey} = $annot_orig{$key};
			#delete($annot_orig{$key});
			$temp_hash_geno{$newkey} = $genos_orig{$key};
			#delete($genos_orig{$key});
			$snps_with_data+=1;
		}
		$key_cnt+=1;
	}
	%annot_orig = %temp_hash_annot;
	%genos_orig = %temp_hash_geno;
	print STDERR "key_cnt =$key_cnt no of SNPs in annotation file.\n";
	print STDERR "$snp_miss_cnt SNPs from the annotation file do not have any genotypes called in the vcf files in \n";
	print STDERR "$snps_with_data SNPs that have both annotation and genotype data\n";
}

sub recode_geno_vcf(){	
my $key;
my %recode_vcfgeno_hash=();
my $alleleA;
my $alleleB;
my $temp;
my $c;
	foreach $key (keys %genos_orig) {
		$alleleA = substr($annot_orig{$key}[2], 0, 1);
		$alleleB = substr($annot_orig{$key}[2], 2, 1);
		$recode_vcfgeno_hash{$alleleA.$alleleA}=0;
		$recode_vcfgeno_hash{$alleleA.$alleleB}=1;
		$recode_vcfgeno_hash{$alleleB.$alleleB}=2;
		for($c=0; $c<$num_files; $c++) {
			if(defined($genos_orig{$key}[$c])) {
				$temp = $genos_orig{$key}[$c];
				$genos_orig{$key}[$c] = $recode_vcfgeno_hash{$genos_orig{$key}[$c]};
				if(!(defined($genos_orig{$key}[$c]))) {
					print STDERR "no call at snp ".$annot_orig{$key}[0]."  alleleA = ".$alleleA."  alleleB = ".$alleleB."  geno = ".$temp."\n";
					$genos_orig{$key}[$c] = -1;
				}
			}
			else{
				$genos_orig{$key}[$c] = -1;
			}
		}
	}
}

sub filter_missing_calls_vcf(){
my $c;
my $key;
my $num_missing;
my $prop_missing;
my $num_orig;
my $num_removed;
	$num_orig = scalar keys %genos_orig;
	print STDERR "Missingness threshold set to ".$missingness.".  Any SNP with more than ".100*$missingness."\% missing calls will be discarded\n";
	print STDERR $num_orig." SNPs with called genotypes prior to filtering for missingness\n";
	foreach $key (keys %genos_orig) {
		$num_missing = 0;
		for($c=0; $c< $num_files; $c++) {
			if($genos_orig{$key}[$c] eq "-1"){
				$num_missing++;
			}
		}
		$prop_missing = $num_missing/$num_files;
		if($prop_missing > $missingness) {
			delete($genos_orig{$key});
			delete($annot_orig{$key});
		}
	}
	$num_removed = $num_orig - scalar keys %genos_orig;
	print STDERR $num_removed." SNPs removed\n";
	print STDERR scalar keys %genos_orig," SNPs remaining\n";
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
	print STDERR $tim;
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
		die "Tar archive reading/writing is not supported. File: $file"; 
	}
	# change mode required for gzip and bzip2 files
	if($file =~ /\.gz$/) {
		# gzip file open
		unless($perliogzip) { die "Perl module PerlIO::gzip not installed. Required for reading/writing gzip compressed files.\nOpen was attempted on: $file" }
		$mode = "$mode:gzip"; 
		if(!defined($copen_gz_warn) && $file =~ /vcf\.gz$/) {
			print STDERR "WARNING: The Perl module PerlIO::gzip cannot read VCF files compressed with the bgzip compression utility and may cause errors.\nIf you have compressed these files with gzip then please ignore this message.\n\n";
			$copen_gz_warn = 1;
		}
	} elsif ($file =~ /\.bz2$/) {
		unless($perlioviabzip2) { die "Perl module PerlIO::via::Bzip2 not installed. Required for reading/writing bzip2 compressed files.\nOpen was attempted on: $file" }
		$mode = "$mode:via(Bzip2)";
	} 
	return( open($_[0], $mode, $file) );
}

