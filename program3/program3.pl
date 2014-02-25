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
    
    my @allMatches = ();
    my @allMatchIndices = ();
    
    for (my $i = 0; $i < scalar(@seqs); $i++) {
        # find the matches in the sequence
        my @matches = $seqs[$i] =~ m/([OUZ]\-+|\-+[OUZ])/g;
        my @matchIndices = ();
        
        # find the indices of the matches in the sequence
        my $currentIndex = 0;
        for (my $j = 0; $j < scalar(@matches); $j++) {
            my $index = index($seqs[$i], $matches[$j], $currentIndex);
            print "$matches[$j] $index\n";
            push (@matchIndices, $index);
            $currentIndex = $index + length($matches[$j]);
        }
        
        push (@allMatches, @matches);
        push (@allMatchIndices, @matchIndices);
    }
    
    my %overlaps = ();

    # go through each collection of matches
    for (my $i = 0; $i < scalar(@allMatches); $i++) {
        my @matches = $allMatches[$i];
        my @matchIndices = $allMatchIndices[$i];
        
        # go through each match in the collection
        for (my $j = 0; $j < scalar(@matches); $j++) {
            my $match = $matches[$j];
            my $matchIndex = $matchIndices[$j];
            #print "$match $matchIndex\n";
            
            # go through the other collections of matches
            for (my $otherI = $i + 1; $otherI < scalar(@allMatches); $otherI++) {
                my @otherMatches = $allMatches[$otherI];
                my @otherIndices = $allMatchIndices[$otherI];
                
                # go through each match in the other collection
                for (my $otherJ = 0; $otherJ < scalar(@otherMatches); $otherJ++) {
                    my $otherMatch = $otherMatches[$otherJ];
                    my $otherMatchIndex = $otherIndices[$otherJ];
                    
                    (my $start, my $end) = &getStringOverlap($match,
                                                             $matchIndex,
                                                             $otherMatch,
                                                             $otherMatchIndex);
                    
                    if ($start != -1 and $end != -1) {
                        my $removed = 0;
                        
                        while ( my ($otherStart, @other) = each(%overlaps) ) {
                            my $otherEnd = $other[0];
                            my $otherCount = $other[1];
                            (my $newStart, my $newEnd) = &getNumberOverlap($start,
                                                                           $end,
                                                                           $otherStart,
                                                                           $otherEnd);
                            if ($newStart != -1 and $newEnd != -1) {
                                delete $overlaps{$otherStart};
                                @{$overlaps{$newStart}} = ($newEnd, $otherCount + 1);  
                                $removed = 1;
                                last;
                            }                                
                        }
                        
                        if ($removed == 0) {
                            @{$overlaps{$start}} = ($end, 1);    
                        }
                    }                    
                }
            }
        }
    }
    
    for my $key ( keys %overlaps ) {
        my $value = $overlaps{$key};
        print "$key => $value\n";
    }
}

# close the log file
close $logFile;

###############################################################################

sub min() {
    my $numA = $_[0];
    my $numB = $_[1];
    
    if ($numA < $numB) {
        return $numA;
    }
    
    return $numB;
}

sub max() {    
    my $numA = $_[0];
    my $numB = $_[1];
    
    if ($numA > $numB) {
        return $numA;
    }
    
    return $numB;
}

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

sub getStringOverlap() {
    my $stringA = $_[0];
    my $startA = $_[1];
    
    my $stringB = $_[2];
    my $startB = $_[3];
    
    my $lengthA = length($stringA);
    my $lengthB = length($stringB);
    
    my $endA = $startA + $lengthA;
    my $endB = $startB + $lengthB;
    
    return &getNumberOverlap($startA, $endA, $startB, $endB); 
}

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