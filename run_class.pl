#!/bin/perl

# The wrapper of running CLASS
# Li Song

use strict ;
use Cwd 'cwd' ;
use Cwd 'abs_path' ;
use File::Basename;

my $CWD = dirname( abs_path( $0 ) ) ;

my $usage = "Usage: perl run_class.pl [options]\n".
				"Options:\n".
				"\t-a alignment_file (REQUIRED): the path to the alignemtn file(in SAM or BAM format)\n".
				"\t-o output_file: the file storing the output of CLASS (default: ./alignment_file_wo_extension.gtf)\n".
				"\t-p number_of_threads: specify the number of worker threads (default:1)\n".
				"\t-F f: do not report the transcripts whose abundance level is lower than f*|most expressed transcript| in a gene\n".
				"\t-l label: add a prefix and a \"_\" to the ids in the GTF file (default: not used)\n".
				"\t-j junction: the path to the splice junction file\n".
				"\t-e evidence: the path to the evidence files\n".
				"\t--var_rd_len: extensive variable read lengths, i.e. reads after trimming (default: no)\n".
				"\t--set-cover: use set cover to build transcripts from splicing graph (default: no)\n".
				"\t--verbose: also output the procedure of CLASS (default: no)\n".
				"\t--wd tempoaray_file_directory: the directory storing the temporary files (default: ./class_tmp)\n".
				"\t--clean: whehter to remove the temporary files in -wd (default: no)\n";


# process the parameters
my $i ;
my $outputFile = "--" ;
my $alignmentFile = "--" ;
my @classARGV ;
my $tmpWd = "./class_tmp" ;
my $cleanTempFiles = 0 ;
my $file ;
my $filePrefix ;
my $alignmentFilePath ;

if ( @ARGV == 0 )
{
	die $usage ;
}

for ( $i = 0 ; $i < @ARGV ; ++$i )
{
	if ( $ARGV[$i] eq "-h" )
	{
		die $usage ;
	}
	elsif ( $ARGV[$i] eq "-o" )
	{
		$outputFile = $ARGV[ $i + 1 ] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "-a" )
	{
		$alignmentFile = $ARGV[$i + 1] ;
		# convert it to full path
		$alignmentFile = abs_path( $alignmentFile ) ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--wd" )
	{
		$tmpWd = $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--clean" )
	{
		$cleanTempFiles = 1 ;
	}
	elsif ( $ARGV[$i] eq "-p" || $ARGV[$i] eq "-F" || $ARGV[$i] eq "-j" || $ARGV[$i] eq "-e" || $ARGV[$i] eq "-l" )
	{
		push @classARGV, $ARGV[$i] ;
		push @classARGV, $ARGV[$i + 1] ;
		++$i ;
	}
	elsif ( $ARGV[$i] eq "--verbose" || $ARGV[$i] eq "--set-cover" || $ARGV[$i] eq "--var_rd_len" ) 
	{
		push @classARGV, $ARGV[$i] ;
	}
	else
	{
		die "Unknown parameter", $ARGV[$i], "\n" ;
	}
}

die "Must specify the alignment file.\n" if ( $alignmentFile eq "--" ) ;
die "The alignment file must be with extension sam or bam" if ( !( $alignmentFile =~ /\.[bs]am$/) ) ;

$file = fileparse( $alignmentFile, qr/$/ ) ;
$filePrefix = fileparse( $alignmentFile, qr/\.[^.]*/) ;
#print $file, " ", $filePrefix, "\n" ;
#exit ; 
if ( $outputFile eq "--" )
{
	$outputFile = "$filePrefix.gtf" ;
}

if ( !(-e $tmpWd) )
{
	`mkdir $tmpWd` ;
}
if ( -e "$tmpWd/$file")
{
	unlink "$tmpWd/$file" ;
}
symlink $alignmentFile, "$tmpWd/$file" ;

print "Generate the depth file.\n" ;
my $skip = 0 ;
if ( -e "$tmpWd/$filePrefix.depth_log" && -e "$tmpWd/$filePrefix.depth" ) # Test whether to skip the step
{
	open FP1, "$tmpWd/$filePrefix.depth_log" ;
	my $line = <FP1> ;
	chomp $line ;
	$skip = 1 if ( $line eq $alignmentFile ) ;
	close FP1 ;
}
if ( $skip == 1 )
{
	print "Found generated depth file. Skip this step.\n" ;
}
else
{
	unlink "$tmpWd/$filePrefix.depth_log" if ( -e "$tmpWd/$filePrefix.depth_log" ) ; # Test whether to skip the step

	print "Running samtools.\n" ;
	#my $pattern = "\"\\t0\$\"" ;
	#print $pattern, "\n" ;
	system( "$CWD/samtools-0.1.19/samtools depth $tmpWd/$file | grep -v -P \"\\t0\$\" > $tmpWd/$filePrefix.depth" ) ;
	#system( "$CWD/samtools-0.1.19/samtools depth $tmpWd/$file > $tmpWd/$filePrefix.depth" ) ;

	open FP1, ">$tmpWd/$filePrefix.depth_log" ;
	print FP1 $alignmentFile ;
	close FP1 ;
}

print "Generate the splice sites file.\n" ;
$skip = 0 ;
if ( -e "$tmpWd/$filePrefix.splice_log" && -e "$tmpWd/$filePrefix.splice" ) # Test whether to skip the step
{
	open FP1, "$tmpWd/$filePrefix.splice_log" ;
	my $line = <FP1> ;
	chomp $line ;
	$skip = 1 if ( $line eq $alignmentFile ) ;
	close FP1 ;
}
if ( $skip == 1 )
{
	print "Found generated splice sites file. Skip this step.\n" ;
}
else
{
	unlink "$tmpWd/$filePrefix.splice_log" if ( -e "$tmpWd/$filePrefix.splice_log" ) ; # Test whether to skip the step

	print "Running junc.\n" ;
	system "$CWD/junc $tmpWd/$file -a > $tmpWd/$filePrefix.splice" ;
	
	open FP1, ">$tmpWd/$filePrefix.splice_log" ;
	print FP1 $alignmentFile ;
	close FP1 ;
}

print "Running CLASS.\n" ;
system "$CWD/class $tmpWd/$filePrefix @classARGV > $outputFile" ;

if ( $cleanTempFiles == 1 )
{
	unlink "$tmpWd/$file", "$tmpWd/$filePrefix.depth", "$tmpWd/$filePrefix.splice", "$tmpWd/$filePrefix.depth_log", "$tmpWd/$filePrefix.splice_log" ;
}

