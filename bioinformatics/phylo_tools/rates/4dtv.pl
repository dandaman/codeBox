#!/usr/bin/perl 
#Created by Cedric Simillion on Fri Jul  1 10:17:12 CEST 2005

=head1 Usage

  4dtv.pl -fasta <file> -pairs <file> [ -method BLAST | CLUSTAL ]

=over

=item -fasta

The name of a fasta file containing CDS (mRNA) sequences.

=item -pairs

The name of a file containing all the pairs for which the K4dtv must be calculated.

=item -method

Optional. Tells which alignment method to use. If you type "BLAST" then each sequence pair
will be aligned using bl2seq (BLAST 2 sequences); if you specify "CLUSTAL" then CLUSTALW
will be used.

The default option is CLUSTAL.

=back

=head1 Description

For each pair of identifiers in the -pairs file, the program aligns their mRNA sequences
per codon using CLUSTAL or bl2seq. Next, the fraction of synonymous transversions at
four-fold degenerate sites (K4dtv) is calculated. The output is a tab delimited file sent
to STDOUT.

=cut

use strict;
use lib qw(/home/lang/scripts/TreePipe/msa);
use bioutils;
use IPC::Open2;

#===============================================================================
# Preparation
#===============================================================================
my ($tmp,$start_path);
my @pairs;
my (%params,%sequences);

%params= @ARGV;

unless ( exists( $params{"-fasta"} ) && exists( $params{"-pairs"} ) )
 {
  system("pod2text $0");
  exit 1;
 } #exists( $params{"-fasta"} ) && exists( $params{"-pairs"} )

exists( $params{"-method"} ) || ( $params{"-method"} = 'CLUSTAL' );

@pairs = &read_pairs( $params{"-pairs"} );
%sequences= &fasta2hash( $params{"-fasta"} );
$tmp = "/tmp/".time.$$;
chomp( $start_path = `pwd` );
mkdir ( $tmp, 0777 ) || die "Could not create temporary directore $tmp : $!\n";
chdir $tmp;

#===============================================================================
# Main loop
#===============================================================================
print STDOUT join("\t", qw(seq1 seq2 L p K) )."\n";
foreach my $pair ( @pairs )
 {

  defined( $sequences{ $$pair[0] } ) || print STDERR "No sequence for $$pair[0].\n";
  defined( $sequences{ $$pair[1] } ) || print STDERR "No sequence for $$pair[1].\n";
  ( defined( $sequences{ $$pair[1] } ) && defined( $sequences{ $$pair[0] } ) ) || next;
  
  foreach my $id ( @$pair )
   {
    $sequences{$id} = uc ( $sequences{$id} );
   } #foreach my $id ( @$pair )

  my $align= ( uc($params{"-method"}) eq 'BLAST' )
   ? &codon_bl2seq( $$pair[0], $sequences{ $$pair[0] }, $$pair[1], $sequences{ $$pair[1] } )
   : &codon_clustalw( $$pair[0], $sequences{ $$pair[0] }, $$pair[1], $sequences{ $$pair[1] } );

  if ( !defined $$align{ $$pair[0] } )
   {
    print STDERR "$$pair[0] and $$pair[1] are not alignable.\n";
    next;
   } #if ( !defined $$align{$ids[0]} )

  &strip_align( $align );
  &fourdtv( $pair, $align, $params{"-runs"} );
 } #foreach my $pair ( @pairs )

chdir $start_path ;
unlink( "$tmp/*" );
rmdir $tmp;

#===============================================================================
# Subroutines
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub read_pairs
# Reads in all pairs from a tab-delimited file. Identical pairs are
# ignored
# Input: a file name.
# Returns: an array or arrays. Each row contains 2 elements
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub read_pairs
 {
  my $file=$_[0];
  my @pairs;
  open( IN, $file ) || die "Could not open file $file: $!\n";
  while ( <IN> )
   {
    chomp;
    my @row;
    @row[0,1]= split( /\t/ );
    ( $row[0] eq $row[1] ) && next;
    push( @pairs, \@row );
   } #while ( <IN> )
  close IN;
  
  return( @pairs );
 } #sub read_pairs

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub codon_bl2seq
# Takes 2 mRNA sequences and translates them in to protein sequences.
# These protein sequences are then aligned to each other using bl2seq.
# The obtained alignment is then back-translated to an alignment of
# the original mRNA sequences. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub codon_bl2seq
 {
  my (@ids, @cds ,@pep, @alignment );
  my (%return);
  ( $ids[0], $cds[0], $ids[1], $cds[1] )=@_;
  @pep = map { &translate($_) } ( @cds );
  
  foreach ( 0, 1 )
   {
    open( OUT, ">seq$_.tfa" ) || die "Could not write sequence file seq$_.tfa: $!\n";
    print OUT ">seq$_\n$pep[$_]\n";
    close OUT;
   } #foreach ( 0, 1 )
  
  open( BLAST, "bl2seq -i seq0.tfa -j seq1.tfa -p blastp -F f -e 0.001 |" );
  while ( <BLAST> ) 
   {

    if ( m/(Query|Sbjct): +(\d+) +(\S+)/ )
     {
      my ($position, $i, $sequence);

      $position = ($2-1)*3;
      $i = ( $1 eq 'Query' ) ? 0  : 1;
      $sequence = $3;

      foreach my $residue ( split( //, $sequence ) )
       {
        if ( $residue ne '-' )
	 {
	  $alignment[$i].= substr( $cds[$i], $position, 3 );
	  $position+= 3;
	 } #if ( $residue ne '-' )
	else 
	 {
	  $alignment[$i].= '---';
	 } #else 
       } #foreach my $residue ( split( //, $sequence ) )

     } #if ( m/Query: +(\d+) +(\S+)/ )

    m/Lambda/ && last;

   } #while ( <BLAST> ) 
  close BLAST;
  
  unlink( "seq0.tfa", "seq1.tfa" );
 
  $return{$ids[0]}=$alignment[0];
  $return{$ids[1]}=$alignment[1]; 

  defined( $return{$ids[0]} ) || print STDERR "$ids[0] and $ids[1] are not alignable.\n";

  return( \%return );
   
 } #sub codon_bl2seq

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub codon_clustalw
# Takes 2 mRNA sequences and translates them in to protein sequences.
# These protein sequences are then aligned to each other using
# CLUSTALW. The obtained alignment is then back-translated to an
# alignment of the original mRNA sequences. 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub codon_clustalw
 {
  my (@ids, @cds ,@pep, @alignment, @position );
  my (%return);
  ( $ids[0], $cds[0], $ids[1], $cds[1] )=@_;
  @pep = map { &translate($_) } ( @cds );
  
  open( OUT, ">seqs.tfa" ) || die "Could not write sequence file seqs.tfa: $!\n";
  foreach ( 0, 1 )
   {
    print OUT ">seq$_\n$pep[$_]\n";
   } #foreach ( 0, 1 )
  close OUT;

  `clustalw seqs.tfa`;

  $position[0]=0;
  $position[1]=0;

  open( IN, "seqs.aln" );
  while( <IN> )
   {
    chomp;
    if ( m/^seq(\d) +(.*)/ )
     {
      my ($i,$sequence)=($1,$2);
      foreach my $residue ( split(//, $sequence ) )
       {
        if ( $residue ne '-' )
	 {
	  while ( &translate( substr( $cds[$i], $position[$i], 3 ) ) eq '*' )
	   {
	    $position[$i]+= 3;
	   } #while ( &translate( substr( $cds[$i], $position[$i], 3 ) ) eq '*' )
	  
	  $alignment[$i].= substr( $cds[$i], $position[$i], 3 );
	  $position[$i]+= 3;
	 } #if ( $residue ne '-' )
	else 
	 {
	  $alignment[$i].= '---';
	 } #else
       } #foreach my $residue ( split(//, $sequence ) )
     } #if ( m/^seq(\d) +(.*)/ )
   } #while( <IN> )
  close IN;

  unlink( "seqs.*" );
  
  $return{$ids[0]}=$alignment[0];
  $return{$ids[1]}=$alignment[1]; 

  defined( $return{$ids[0]} ) || print STDERR "$ids[0] and $ids[1] are not alignable.\n";

  return( \%return );

 } #sub codon_clustalw

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub strip_align
# Calls the strip_align2 program to strip gaps from an alignment.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub strip_align
 {
  my ($align, $read, $write, $id, $process_id);
  $align= $_[0];
#  $process_id = open2( $read, $write, "~/scripts/TreePipe/msa/stripalign3.pl" );
  $process_id = open2( $read, $write, "~/scripts/TreePipe/msa/stripalign2.pl 0 100 y" );

#   open( BUG, ">test.tfa" ); #debug
  foreach $id (keys %$align)
   {
#     print BUG ">$id\n$$align{$id}\n"; #debug
    print $write ">$id\n$$align{$id}\n";
    delete $$align{$id};
   } #foreach my $id (keys %$align)
  close $write;
#   close BUG; #debug
  
  while ( <$read> )
   {
    chomp;
    if ( m/^>(.*)/ )
     {
      $id = $1;
     } #if ( m/^>(.*)/ )
    else 
     {
      $$align{$id}.=$_;
     } #else 
   } #while ( <$read> )
  close $read; 
  
  waitpid( $process_id, 0 );
   
 } #sub strip_align

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub fourdtv
# Takes an alignment of 2 CDS's and computes the K4dtv value for it.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub fourdtv
 {
  my ($align,$identities,$differences,$p,$K,$L);
  my (@ids,@sequences);
  my (%is_4d_codon, %group);
  
  %is_4d_codon = ( 'TC' => -1, 'CT' => -1, 'CC' => -1, 'CG' => -1, 'AC' => -1, 'GT' => -1, 'GC' => -1, 'GG' => -1 );
  %group = ( 'C' => 'Y', 'T' => 'Y', 'A' => 'R', 'G' => 'R' );
  
  @ids = @{$_[0]};
  $align= $_[1];
  @sequences = values %$align;
  
  $identities = 0;
  $differences = 0;
  
  for ( my $i=0; $i < length( $sequences[0] ); $i+= 3 )
   {
    my @first_two = map { substr( uc($_), $i, 2 ) } ( @sequences );
    my @third_position = map { substr( uc($_), $i+2, 1 ) } ( @sequences );
    if ( $is_4d_codon{ $first_two[0] } && ( $first_two[0] eq $first_two[1] ) )
     {
      ( defined( $group{ $third_position[0] } ) && defined( $group{ $third_position[1] } ) ) || next;
      ( $group{ $third_position[0] } ne $group{ $third_position[1] } ) ? $differences++ : $identities++;
     } #if ( $is_4d_codon{ $first_two[0] } && ( $first_two[0] eq $first_two[1] ) )
   } #for ( my $i=0; $i < length( $sequences[0] ); $i+= 3 )
  
  $L = $differences + $identities;
  
  if ( $L > 0 )
   {
    $p = $differences/($L);
    $K = ($p < 0.5) ? -1/2*log( 1-2*$p ) : "---";
   } #if ( $L > 0 )
  else 
   {
    $p = '---';
    $K = '---';
   } #else 

  print STDOUT join("\t", @ids, $L, $p, $K)."\n";
 } #sub fourdtv
