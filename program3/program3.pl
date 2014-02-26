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
my $outputDir = &getOutputDirectory($inputDir);

# open the log file
my $logFileBase = $0;
$logFileBase =~ s/.pl$//; 
my $logFileName = &getLogFileName($logFileBase, 0);
open my $logFile, ">", $logFileName or die "Problem opening file $logFileName";

# go through each file in input directory
my @files = glob "$inputDir/*.cut";

my @headers = ();
my @seqs = ();

for my $filename (@files) {
    @headers = ();
    @seqs = ();
    
    print "$filename\n";
    
    open (IN, $filename) or die "Unable to open: $filename";
    
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
    
    my %overlaps = ();
    
    for (my $i = 0; $i < scalar(@seqs); $i++) {
        # find the matches in the sequence
        my @matches = $seqs[$i] =~ m/([OUZ]\-+|\-+[OUZ])/g;
        my @matchIndices = ();
        
        # find the indices of the matches in the sequence
        my $currentIndex = 0;
        for (my $j = 0; $j < scalar(@matches); $j++) {
            my $match = $matches[$j];
            my $start = index($seqs[$i], $match, $currentIndex);
            my $end = $start + length($match);
            
            (my $letter = $match) =~ s/-//g;
            
            my $existingOverlapModified = 0;
                        
            for my $otherStart ( keys %overlaps ) {                            
                (my $otherEnd, my $otherCount, my $ouzCode) = @{$overlaps{$otherStart}};
                
                (my $newStart, my $newEnd) = &getNumberOverlap($start,
                                                               $end,
                                                               $otherStart,
                                                               $otherEnd);
                if ($newStart != -1 and $newEnd != -1) {
                    my $newOuz = &buildOuzCode($ouzCode, $letter, $i);;
                                        
                    delete $overlaps{$otherStart};
                    @{$overlaps{$newStart}} = ($newEnd, $otherCount + 1, $newOuz);
                    $existingOverlapModified = 1;
                    last;
                }                                
            }
            
            if (!$existingOverlapModified) {
                my $ouzCode = &buildOuzCode("...", $letter, $i);
                
                @{$overlaps{$start}} = ($end, 1, $ouzCode);    
            }
            
            $currentIndex = $start + length($matches[$j]);
        }
    }
    
    for my $start ( sort { $a <=> $b} keys %overlaps ) {
        my @values = @{$overlaps{$start}};
        my $end = $values[0];
        my $count = $values[1];
        my $ouzCode = $values[2];
        
        if ($count > 1) {
            print "$start => $end, $count, $ouzCode\n";
            
            &writeToLog($headers[0], $start, $end, $ouzCode, "", "");
        }
    }
}

# close the log file
close $logFile;

###############################################################################

sub writeToLog() {
    my $orthId = $_[0];
    my $start = $_[1];
    my $stop = $_[2];
    my $ouzCode = $_[3];
    my $revisionFlag = $_[4];
    my $comment = "";
    
    if (@_ == 6) {
        $comment = $_[5];    
    }
    
    my $line = "$orthId\t$start\t$stop\t$ouzCode\t$revisionFlag\t$comment";

    print{$logFile} "$line\n";
}

# Builds a new OUZ string.
#
#   oldOuz - the OUZ code to be modified
#   letter - the letter that is being added to the code
#   seqIndex - the index of the gene sequence
#
sub buildOuzCode() {
    my $oldOuzCode = $_[0];
    my $letter = $_[1];
    my $seqIndex = $_[2];
    
    if ($seqIndex == 0) {
        return $letter . $oldOuzCode;
    }
        
    if ($seqIndex == 3) {
        return $oldOuzCode . $letter;
    }
        
    my $firstPeriod = index($oldOuzCode, "\.");
    my $secondPeriod = index($oldOuzCode, "\.", $firstPeriod + 1);
    
    my $pre;
    my $post;
    
    if ($seqIndex == 1) {
        $pre = substr($oldOuzCode, 0, $firstPeriod + 1);
        $post = substr($oldOuzCode, $firstPeriod + 1);
    }
    
    if ($seqIndex == 2) {
        $pre = substr($oldOuzCode, 0, $secondPeriod + 1);
        $post = substr($oldOuzCode, $secondPeriod + 1);
    }
    
    return $pre . $letter . $post;
}

# Returns the minimum of two numbers.
# 
#   numA - the first number
#   numB - the second number
# 
sub min() {
    my $numA = $_[0];
    my $numB = $_[1];
    
    if ($numA < $numB) {
        return $numA;
    }
    
    return $numB;
}

# Returns the maximum of two numbers.
# 
#   numA - the first number
#   numB - the second number
# 
sub max() {    
    my $numA = $_[0];
    my $numB = $_[1];
    
    if ($numA > $numB) {
        return $numA;
    }
    
    return $numB;
}

# Find the overlap between two intervals if it exists, otherwise returns
# (-1, -1).
#
#   startA - the start of the first interval
#   endA - the end of the first interval
#   startB - the start of the second interval
#   endB - the end of the second interval
#
sub getNumberOverlap() {
    my $startA = $_[0];
    my $endA = $_[1];
    
    my $startB = $_[2];
    my $endB = $_[3];
    
    # A ends before B starts
    if ($endA < $startB) {
        return (-1, -1); # no overlap
    }
    
    # B ends before A starts
    if ($endB < $startA) {
        return (-1, -1); # no overlap
    }
    
    my $start = &min($startA, $startB);
    my $end = &max($endA, $endB);
    
    return ($start, $end);  
}

# Gets the output directory from the command line argument if it was specified
# or makes the output directory based on the input directory.  Creates the
# directory if it does not exist.
# 
sub getOutputDirectory() {
    my $inputDirectory = $_[0];
    my $outputDir = "$inputDirectory-mod";
    
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

# Gets the input directory from the command line argument and determines
# whether or not the input directory exists.
#
sub getInputDirectory() {
    my $inputDir = $ARGV[ 0 ];
    
    # make sure the specified input directory exists
    if (!-d $inputDir) {
        die "$inputDir is not a directory.";
    }
    
    return $inputDir;
}

# Recursively constructs the file name.
#
#   fileBase - the base of the log file name
#   fileNumber - the current number to attempt to make a file with
#
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
#
sub checkUsage() {
   if ( @ARGV == 0 || $ARGV[0] eq "-h" ) {
      print STDERR "$usageMsg\n";
      exit;
   }
}