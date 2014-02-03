#!/usr/bin/perl -w
#
# Carmen St. Jean (crr8@unh.edu)
# MCBS 913, Spring 2014
# Program 2
#
# February 3, 2014
#

use warnings;
use strict;

package main;

my $aminoStop = "*";
my $nucleoUnknown = "N";

my %aminoAcids = (
    'UCA' => 'S', # Serine
    'UCC' => 'S', # Serine
    'UCG' => 'S', # Serine
    'UCU' => 'S', # Serine
    'UCN' => 'S', # Serine
    
    'UUC' => 'F', # Phenylalanine
    'UUU' => 'F', # Phenylalanine
    'UUA' => 'L', # Leucine
    'UUG' => 'L', # Leucine
    
    'UAC' => 'Y', # Tyrosine
    'UAU' => 'Y', # Tyrosine
    'UAA' => $aminoStop, # Stop
    'UAG' => $aminoStop, # Stop
    
    'UGC' => 'C', # Cysteine
    'UGU' => 'C', # Cysteine
    'UGA' => $aminoStop, # Stop
    'UGG' => 'W', # Tryptophan
    
    'CUA' => 'L', # Leucine
    'CUC' => 'L', # Leucine
    'CUG' => 'L', # Leucine
    'CUU' => 'L', # Leucine
    'CUN' => 'L', # Leucine
    
    'CAU' => 'H', # Histidine
    'CAA' => 'Q', # Glutamine
    'CAG' => 'Q', # Glutamine
    'CAC' => 'H', # Histidine
    
    'CGA' => 'R', # Arginine
    'CGC' => 'R', # Arginine
    'CGG' => 'R', # Arginine
    'CGU' => 'R', # Arginine
    'CGN' => 'R', # Arginine
    
    'AUA' => 'I', # Isoleucine
    'AUC' => 'I', # Isoleucine
    'AUU' => 'I', # Isoleucine
    'AUG' => 'M', # Methionine
    
    'ACA' => 'T', # Threonine
    'ACC' => 'T', # Threonine
    'ACG' => 'T', # Threonine
    'ACU' => 'T', # Threonine
    'ACN' => 'T', # Threonine
    
    'AAC' => 'N', # Asparagine
    'AAU' => 'N', # Asparagine
    'AAA' => 'K', # Lysine
    'AAG' => 'K', # Lysine
    
    'AGC' => 'S', # Serine
    'AGU' => 'S', # Serine
    'AGA' => 'R', # Arginine
    'AGG' => 'R', # Arginine
    
    'CCA' => 'P', # Proline
    'CCC' => 'P', # Proline
    'CCG' => 'P', # Proline
    'CCU' => 'P', # Proline 
    'CCN' => 'P', # Proline     
    
    'GUA' => 'V', # Valine
    'GUC' => 'V', # Valine
    'GUG' => 'V', # Valine
    'GUU' => 'V', # Valine
    'GUN' => 'V', # Valine
    
    'GCA' => 'A', # Alanine
    'GCC' => 'A', # Alanine
    'GCG' => 'A', # Alanine
    'GCU' => 'A', # Alanine
    'GCN' => 'A', # Alanine
    
    'GAC' => 'D', # Aspartic Acid
    'GAU' => 'D', # Aspartic Acid
    'GAA' => 'E', # Glutamic Acid
    'GAG' => 'E', # Glutamic Acid
    
    'GGA' => 'G', # Glycine
    'GGC' => 'G', # Glycine
    'GGG' => 'G', # Glycine
    'GGU' => 'G', # Glycine
    'GGN' => 'G'  # Glycine
);

my $usageMsg = q(   Usage: program2 fastafile

          Extract each sequence from a fastafile into a single string.
          <do something to the sequence -- this one computes its length
          and adds it after the sequence name on the header>

          Output is the revised header and sequence data
          Output sent to standard output. );
          
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
&checkUsage();

my $seqFile = $ARGV[ 0 ];

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
   my $rna = &dnaToRna($seq);
   my $aminoAcids = &rnaToAminoAcids($rna);
   
   print "$header\n$aminoAcids\n\n";
   #--------------------------------------------------------
   $header = $inLine;    # last line read is either next header or null
}

#my $dna = "ATCATCATCATCGGGTGA";
#my $rna = &dnaToRna($dna);
#my $rev = &reverseComplement($dna);
#my $aa = &rnaToAminoAcids($rna);

#print "rev: $rev\n";
#print "rna: $rna\n";
#print "dna: $dna\n";
#print "aa: $aa\n";

# Translates an RNA sequence to an amino acid sequence.
sub rnaToAminoAcids() {
    my $rna = $_[0];
    my $protein = "";
    
    for (my $i = 0; ($i + 3) <= length($rna); $i += 3) {
        my $codon = substr($rna, $i, 3);
        my $aminoAcid = &codonToAminoAcid($codon);
        
        if ($aminoAcid) {
            $protein .= $aminoAcid;
        } else {
            print "\n***\n$codon\n***\n";
            
            #$protein .= "-";
        }
    }
    
    return $protein;
}



# Takes an RNA codon and returns an amino acid.
sub codonToAminoAcid() {
    my $codon = $_[0];
    
    return $aminoAcids{$codon};
}

# Transcribes from DNA to RNA.
sub dnaToRna() {
    our $rna = $_[0];
    
    $rna =~ tr/T/U/;
    
    return $rna;
}

# Creates the reverse complement of the DNA.
sub reverseComplement() {
    our $dna = $_[0];
    
    our $revCmp = reverse $dna;
    $revCmp =~ tr/ACTG/TGAC/;
    
    return $revCmp;
}

# Checks usage.
sub checkUsage() {
   if ( @ARGV == 0 || $ARGV[0] eq "-h" ) {
      print STDERR "$usageMsg\n";
      exit;
   }
}