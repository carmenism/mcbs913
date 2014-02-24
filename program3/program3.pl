#!/usr/bin/perl -w
#
# Carmen St. Jean (crr8@unh.edu)
# MCBS 913, Spring 2014
# Program 3
#
# February 24, 2014
#

use warnings;
use strict;

my $usageMsg = q(   Usage: program3.pl fastafile

          Extract each sequence from a fastafile into a single string.
          Transcribes the DNA to RNA, then translates RNA codons to
          amino acids.

          Output is the revised header and protein sequences.
          Output sent to standard output. );

###############################################################################

# check usage and open the input file
&checkUsage();
my $inputDir = &getInputDirectory();
my $outputDir = &getOutputDirectory();

# open the log file
my $logFileBase = $0;
$logFileBase =~ s/.pl$//; 
my $logFileName = &getLogFileName($logFileBase, 0);
open my $logFile, ">", $logFileName or die "Problem opening file $logFileName";

# go through each file in input directory
my @files = glob "$inputDir/*";

my @headers = ();
my @seqs = ();

for my $filename (@files) {
    @headers = ();
    @seqs = ();
    
    print "$filename\n";
    
    open (IN, $filename) or die "Unable to open: " . $filename;
    
    # first line better be a sequence header
    my $header = <IN>;
    if ( substr( $header, 0, 1 ) ne '>' ) {
        print "********* ERROR ********* Expecting header, got:\n $header";
        print " is this a fasta file??? ";
        &checkUsage();
        exit;
    }
    
    while ($header) {
        my $seq = ""; 
        my $inLine = <IN>;
    
        # read in all input lines of bases
        while ( $inLine && substr( $inLine, 0, 1 ) ne '>' ) {
            chomp( $inLine );     # remove line feed
            $seq = $seq . $inLine;
            $inLine = <IN>;
        }
    
        chomp( $header );  # remove line feed
        $header .= " ";    # make sure there is at least one space after seq id
        my $seqId = substr( $header, 1, index( $header, " " ) - 1 );
        
        push (@headers, $header);
        push (@seqs, $seq);
        
        $header = $inLine;
    }
    
    # DO WORK HERE
    
    exit();
}

# close the log file
close $logFile;

###############################################################################

sub getOutputDirectory() {
    my $outputDir = "$inputDir-mod";
    
    # check if the output directory parameter exists
    if (@ARGV == 2) {
        $outputDir = $ARGV[1];
    }
    
    # make the output directory if it does not exist
    if (!-d $outputDir) {
        mkdir($outputDir, 0700)
        or die "Problem creating output directory $outputDir.";
    }
}

sub getInputDirectory() {
    my $inputDir = $ARGV[ 0 ];
    
    # make sure the specified input directory exists
    if (!-d $inputDir) {
        die "$inputDir is not a directory.";
    }
    
    return $inputDir;
}

# Recursively constructs the file name.
sub getLogFileName() {
    my $fileBase = $_[0];
    my $fileNumber = $_[1];    
    my $name;
    
    if ($fileNumber == 0) {
        $name = "$fileBase.log";
    } else {
        $name = "$fileBase-$fileNumber.log";
    }
    
    if (!-e $name) {
        return $name;
    }
    
    return &getLogFileName($fileBase, ($fileNumber + 1));
}

# Checks usage.
sub checkUsage() {
   if ( @ARGV == 0 || $ARGV[0] eq "-h" ) {
      print STDERR "$usageMsg\n";
      exit;
   }
}