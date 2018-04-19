#!/usr/bin/perl 
#Created by Cedric Simillion on Thu Oct  3 15:35:56 MEST 2002

=head1 Usage
 
 > stripalign2.pl < my_alignment.tfa <gaps_tolerated> <p> [codon]

Reads an alignment (fasta format) from STDIN. Strips every position in the alignment that
is considered as a gap position + all positions from both ends of the gap until the next
"anchor" positions on both sides.

A column is considered as a gap position if more than a fraction <gaps_tolerated> (a value
between 0 and 1) of the positions in this column are gaps. An anchor position is a column
where the residues are more or less conserved. This is determined as follows: for every
pair of residues in the column, the BLOSUM62 value is retrieved. Next, the <p>th
percentile for all these values is calculated. If this percentile value is bigger or equal
than zero, the column is considered as an anchor position. In other words, the <p>% lowest
BLOSUM62 values must be non-negative values, which is in accordance to the BLOSUM62
matrix, where a non-negative value for a given amino acid pair indicates that this
substitution is likely to happen.

When working with mRNA sequences, the program can edit the alignment per codon by
specifying the parameter (just type something, e.g. 'y' as parameter). Anchor positions
will then be defined as positions were all the codons code for the same amino acid. It is
then assumed that the protein translation is in reading frame +1.

The stripped alignment is returned to STDOUT in fasta format

=cut

use strict;
use lib qw(/home/pgsb/daniel.lang/git/codebox/bioinformatics/phylo_tools/rates);
use bioutils;
use Statistics::Descriptive;

#===============================================================================
# Initialisation
#===============================================================================
my ($per_codon,$c,$i,$last_anchor,$stripping,$max_gap_portion,$percentile);
my (@align,@ids);
my (%align,%blosum);
($max_gap_portion,$percentile,$per_codon)=@ARGV;
($per_codon) ? ($per_codon=3) : ($per_codon=1);
%align=&input2hash;
%blosum=&read_blosum;

#===============================================================================
# Main
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Construct a data structure for editing the alignment
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$c=0;
foreach (keys %align)
 {
  my ($n);
  $ids[$c]=$_;
  $n=0;
  for ($i=0;$i<length($align{$_});$i+=$per_codon)
   {
    $align[$n][$c]=substr($align{$_},$i,$per_codon);
    $n++;
   } #for ($i=0;$i<length($align{$_};$i+=$per_codon)
  $c++; 
 } #foreach (keys %align)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Strip the alignment
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$last_anchor=-1;
for ( $i=0;$i<scalar(@align);$i++ )
 {
  &anchor( @{$align[$i]} ) && ( $last_anchor=$i );
  if ( &gap( @{$align[$i]} ) || ( $last_anchor == -1 ) )
   {
#    while ( !&anchor( @{$align[$last_anchor]} ) && ( $last_anchor > -1 ) ) {$last_anchor--}
    $last_anchor++;
    do
     {
      splice(@align,$last_anchor,1);
     } until ( ($last_anchor == scalar(@align)) || &anchor(@{$align[$last_anchor]}) );
    $i=$last_anchor;
   } #if (&gap(@{$align[$i]))
 } #for ($i=0;$i<scalar(@align);$i++)
( $last_anchor < $#align ) && splice (@align,$last_anchor+1); 
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Store the edited alignment again in the original hash structure
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%align=();
foreach (@align)
 {
  for ($i=0;$i<scalar(@{$_});$i++)
   {
    $align{$ids[$i]}.=$$_[$i];
   } #for ($i=0;$i<scalar(@{$_});$i++)
 } #foreach (@align)

foreach (keys %align)
 {
  print STDOUT ">$_\n$align{$_}\n";
 } #foreach (keys %align)

#===============================================================================
# Subroutines
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub fasta2hash
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub input2hash
 {
  my ($key);
  my (%fasta_hash);
  while (<STDIN>)
   {
    chomp;
    if (/^>(\S+)/)
     {
      $key=$1;
     } #if (/^>(\w)$/)
    else
     {
      $key || die "Input is not in fasta format!";
      s/\s+//g;
      s/[^ARNDCQEGHILKMFPSTWYV*-]/X/gi;
      $fasta_hash{$key}.=$_;
     } #else
   } #while (<STDIN>)
  return (%fasta_hash);
 } #sub input2hash ( $ )

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub gap
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub gap
 {
  my ($test,$gaps);
  $gaps=0;
  $test= '-' x $per_codon;
  foreach (@_)
   {
    ($_ eq $test) && ($gaps++);
   } #foreach (@_)
  ($gaps/scalar(@_)>$max_gap_portion) ? (return -1) : (return 0); 
 } #sub gap

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub anchor
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
sub anchor
 {
  my (@row)=@_;
  my ($return,$test,$gap);
  $return=0;
  $gap='-' x $per_codon;
  unless (&gap(@row))
   {
    my ($i,$blosums);
    $blosums=Statistics::Descriptive::Full->new();
    if ($per_codon==3)
     {
      @row=map  {$_=&translate($_)}(@row);
     } #if ($per_codon==3)

    if ( scalar( @row  ) > 2 )
     {
      foreach $i (0..$#row-1)
       {
	my ($j);
	($row[$i] eq $gap) && (next);
	foreach $j ($i+1..$#row)
	 {
	  ($row[$j] eq $gap) && (next);
          $blosums->add_data($blosum{$row[$i]}{$row[$j]});
	 } #foreach $j ($i+1..$#row)
       } #foreach $i (0..$#row)
      ($blosums->percentile($percentile)>=0) && ($return=-1);
     } #if ( scalar( @row  ) > 2 )
    else 
     {
      ( $blosum{$row[0]}{$row[1]} >= 0 ) && ($return=-1);
     } #else 
   } # unless (&gap(@row))  
  return ($return); 
 } #sub anchor

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub read_blosum
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub read_blosum
 {
  my $blosum_file="/home/pgsb/daniel.lang/git/codebox/bioinformatics/phylo_tools/rates/blosum62.txt";
  my (@residus);
  my (%matrix);
  open (IN,$blosum_file) || die "Could not open $blosum_file: $!\n";
  (undef,@residus)=split(/[\t\n]/,<IN>);
  while (<IN>)
   {
    chomp;
    my ($this,@row)=split(/\t/);
    @{$matrix{$this}}{@residus}=@row;
   } #while (<IN>)
  close IN; 
  return (%matrix) 
 } #sub read_blosum
