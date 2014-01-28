#!/usr/bin/perl -w
#
# fastaParse fastafile
#
#    
#   rdb
#   01/28/10
#
my $usageMsg = q(   Usage: fastaparse fastafile

          Extract each sequence from a fastafile into a single string.
          <do something to the sequence -- this one computes its length
          and adds it after the sequence name on the header>

          Output is the revised header and sequence data
          Output sent to standard output. );

use warnings;
use strict;

#++++++++++++++++++++++ set pattern here ++++++++++++++++++++++++++++
#
#my $pattern = "(AT*).*(TAG|ATG)";

print "Enter a pattern> \n";

my $pattern = <STDIN>;
chomp( $pattern );

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
&checkUsage();              # comment this line if assigning file name above

my $seqFile = $ARGV[ 0 ];   # comment this line if assigning file name above

open ( IN, $seqFile )  or die "Unable to open: ".$seqFile ;

# first line better be a sequence header
my $header = <IN>;
if ( substr( $header, 0, 1 ) ne '>' ) {
   print "********* ERROR ********* Expecting header, got:\n $header";
   print " is this a fasta file??? ";
   &checkUsage();
   exit;
}

while ( $header ) {
   my $seq = ""; 
   my $inLine = <IN>;

   # read in all input lines of bases
   while ( $inLine && substr( $inLine, 0, 1 ) ne '>' ) {
      chomp( $inLine );     # remove line feed
      $seq = $seq . $inLine;
      $inLine = <IN>;
   }

   # -----------------------------------------------------
   #   Replace the lines below with the sequence specific
   #  processing you want to do.
   #
   my $basesCount = length( $seq );

   chomp( $header );  # remove line feed
   $header .= " ";    # make sure there is at least one space after seq id
   #
   # Extract the sequence id field: everything up to the first space.
   #   don't include the '>'; subtract 1 from position of space since
   #   the index includes the '>', but the substring doesn't
   my $seqId = substr( $header, 1, index( $header, " " ) - 1 );

   #print "Sequence: $seqId has a length of $basesCount\n";
   
   my @patternMatches = $seq =~ m/$pattern/g;
   my $patternCount = @patternMatches;
   
   print "Sequence: $seqId, Pattern: $pattern, Number of Matches: $patternCount\n";
   #--------------------------------------------------------
   $header = $inLine;    # last line read is either next header or null
}
#+++++++++++++++++++++++++++++++++++++++++++++++
#                checkUsage
#
sub checkUsage() {
   if ( @ARGV == 0 || $ARGV[0] eq "-h" ) {
      print STDERR "$usageMsg\n";
      exit;
   }
}
