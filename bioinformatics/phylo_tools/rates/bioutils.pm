package bioutils;

use strict;

=head1 Description

  Contains some usefull functions such as

=head2 reverse(seq)

  seq is a string, this function returns another string with all characters in reverse order

=head2 complement(seq)

  seq is a flat nucleic acid sequence stored in a string, 
  this function returns the complement of the sequence (in uppercase)

=head2 translate(seq)

  seq is a flat nucleic acid sequence stored in a string, 
  this function returns the proteic translation of this sequence (in uppercase)

=head2 fasta2flat(file) fasta2flat(file, AC)

  file is the path of a FASTA formated file.
  This function returns the first flat sequence of this file or the one with the given 
  accession number AC (an accession number is a valid identifier WITHOUT internal blanks)

=head2 flat2fasta(file, AC, seq, line_length)

  This function APPENDS a sequence to a FASTA formatted file (file), 
  creates a header line with AC and formats the sequence (seq) with line_length375
  characters per lines

=head2 fasta_cut(file, begin, end)

  This function reads a ONE ENTRY FASTA file and return the sub sequence from begin to end
  
=head2 fasta_length(file)

  This function reads a ONE ENTRY FASTA file and return the length of the sequence

=head2 flat2gb(file, LOCUS, seq)

  This function CREATE a GenBank formated file containing the sequence
  given as parameter, and labelled with LOCUS.

=head2 fasta2hash(file)

  Reads in the given (multiple) fasta file and returns a hash in which the keys are the 
  identifiers of the identifiers of the sequences (without the '>') and the values are 
  the sequences themselves. Created by cesim 14/01/2002

=cut
  
BEGIN {
  use Exporter ();
  use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
  
  # set version
  $VERSION = 1.00;
  
  @ISA = qw(Exporter);
  @EXPORT = qw(&reverse &complement &translate &fasta2flat 
  	       &flat2fasta &gc_percent &flat2gb &fasta_cut &fasta_length &fasta2hash);
  %EXPORT_TAGS = ();
  
  @EXPORT_OK = qw();
}

use vars @EXPORT_OK;


#------------------------------------------------------------
# reverse a flat sequence
sub reverse ( $ )
  {
    my ($seq) = @_;

    return join('', reverse (split '', $seq));
  }

#------------------------------------------------------------
# complement a flat sequence
sub complement( $ )
  {
    my ($seq) = @_;
    
    $seq = uc($seq);
    $seq =~ tr/('A','T','G','C')/('T','A','C','G')/;
    return $seq;
  }

#------------------------------------------------------------
# return the GC% of a sequence (float between 0 and 100)

sub gc_percent( $ )
{
  my($s) = @_;
  (length($s)>0) || die "try to compute GC% on empty sequence";
  return (($s =~ s/[gcGC]//g) * 100) / length($_[0]);
}  

#------------------------------------------------------------
# translate a flat sequence

sub translate ( $ )
{
  my($seq) = @_;
  my($i,$len,$output) = (0,0,'');
  my($codon) = "";

  for($len=length($seq),$seq = uc($seq),$i=0; $i<($len-2) ; $i+=3) {
    $codon = substr($seq,$i,3);

    # would this be easier with a hash system (?) EB

    if   ($codon =~ /^TC/)     {$output .= 'S'; }       # Serine
    elsif($codon =~ /^TT[TC]/) {$output .= 'F'; }       # Phenylalanine
    elsif($codon =~ /^TT[AG]/) {$output .= 'L'; }       # Leucine
    elsif($codon =~ /^TA[TC]/) {$output .= 'Y'; }       # Tyrosine
    elsif($codon =~ /^TA[AG]/) {$output .= '*'; }       # Stop
    elsif($codon =~ /^TG[TC]/) {$output .= 'C'; }       # Cysteine
    elsif($codon =~ /^TGA/)    {$output .= '*'; }       # Stop
    elsif($codon =~ /^TGG/)    {$output .= 'W'; }       # Tryptophan
    elsif($codon =~ /^CT/)     {$output .= 'L'; }       # Leucine
    elsif($codon =~ /^CC/)     {$output .= 'P'; }       # Proline
    elsif($codon =~ /^CA[TC]/) {$output .= 'H'; }       # Histidine
    elsif($codon =~ /^CA[AG]/) {$output .= 'Q'; }       # Glutamine
    elsif($codon =~ /^CG/)     {$output .= 'R'; }       # Arginine
    elsif($codon =~ /^AT[TCA]/){$output .= 'I'; }       # Isoleucine
    elsif($codon =~ /^ATG/)    {$output .= 'M'; }       # Methionine
    elsif($codon =~ /^AC/)     {$output .= 'T'; }       # Threonine
    elsif($codon =~ /^AA[TC]/) {$output .= 'N'; }       # Asparagine
    elsif($codon =~ /^AA[AG]/) {$output .= 'K'; }       # Lysine
    elsif($codon =~ /^AG[TC]/) {$output .= 'S'; }       # Serine
    elsif($codon =~ /^AG[AG]/) {$output .= 'R'; }       # Arginine
    elsif($codon =~ /^GT/)     {$output .= 'V'; }       # Valine
    elsif($codon =~ /^GC/)     {$output .= 'A'; }       # Alanine
    elsif($codon =~ /^GA[TC]/) {$output .= 'D'; }       # Aspartic Acid
    elsif($codon =~ /^GA[AG]/) {$output .= 'E'; }       # Glutamic Acid
    elsif($codon =~ /^GG/)     {$output .= 'G'; }       # Glycine
    else {$output .= 'X'; }                        # Unknown Codon
  }

  return $output;
}


#------------------------------------------------------------
# reas a FASTA file and return the first flat sequence if no AC is given
# or the one that have the given AC (empty if none match)

sub fasta2flat
  {
    my($fastaF, $AC);
    my ($ac, $seq, $inseq) = ('', '', 0);

    if (scalar @_ > 0)
      {
	$fastaF = $_[0];
	if (scalar @_ == 1) {$AC = '';} 
	else {$AC = $_[1];}
      }
    else { die "usage: fasta2flat file_name [AC]";}

    open(FF,$fastaF) || die "Can't open $fastaF";
    while (<FF>)
      {
	$_ =~ s/[\n\r\f]//g;

	if (/^>/)
	  {
	    if ($inseq == 1) {last;} # end of the sequence
	    ($ac) = (/^>([^\s]+)/); # retrieve current accession number
	    if ($AC eq '' || $AC eq $ac) # it's the sequence we are looking for
	      {
		$inseq = 1;
		next;
	      }  
	  }

	if ($inseq == 1)
	  {
	    $seq .= $_;
	    next;
	  }
      }
    close(FF);
    return $seq;
  }

#------------------------------------------------------------
# append a sequence to a FASTA file
#
sub flat2fasta ( $ $ $ $ )
  {
    my($fname, $AC, $seq, $line_length) = @_;
    my($i);
    
    open(FF,">>$fname") || die "can't write into $fname";
    print FF ">$AC\n";
    for ($i=0; $i < length($seq); $i += $line_length)
      {
	printf FF "%s\n", substr($seq,$i,$line_length);
      }
    close(FF);
  }


#------------------------------------------------------------
# extract a sub sequence from a fasta file
#

sub fasta_cut ( $ $ $ )
  {
    my($fname, $begin, $end) = @_;
    my($seq, $i, $copy, $stop, $L, $l, $ac, $b, $e);
    
    ($begin =~ /\d+/) || die "fasta_cut: begin is not a number\n";
    ($end =~ /\d+/) || die "fasta_cut: end is not a number\n";
    ($end >= $begin) || die "fasta_cut: end is greater than begin\n";

    $seq = '';
    $copy=$stop=0;
    $L=0;

    open(FILE_IN,$fname) || die "fasta_cut: can't open file $fname";
    $_ = <FILE_IN>; # header line
    (($ac) = (/^>([^\s]+)/)) || die "fasta_cut: $fname is not a FASTA file (no header)";

    while (<FILE_IN>) # read sequence
      {
	chomp;
	$_ =~ s/\s//g;
	($_ !~ /^>/) || die "fasta_cut: reach end of sequence $fname";
	$l = length($_);
	if ($begin > $L && $begin <= ($L + $l))
	  {
	    $b = $begin - $L - 1;
	    $copy=1;
	  }
	else
	  {
	    $b = 0;
	  }
	if ($end > $L && $end <= ($L + $l))
	  {
	    $e = $end - $L - 1;
	    $stop=1;
	  }
	else
	  {
	    $e = $l - 1;
	  }
	if ($copy==1)
	  {
	    $seq .= substr($_, $b, $e - $b + 1);
	    if ($stop==1) {last;};
	  }

	$L += $l;
      }

    close(FILE_IN);
    return($seq);
  }

#------------------------------------------------------------
# return sequence length
#

sub fasta_length ( $ )
  {
    my ($n);
    open(FIN, $_[0]) || die "Can't open $_[0]";
    $n=-1;
    while(<FIN>)
    {
      if ($n == -1)
      {
	($_ =~ /^>/) || die "$_[0] is not a FASTA file";
	$n =0; next;
      }
      if (/^>/) { last; }
      s/\s+//g;
      $n += length($_);
    }
    close(FIN);

    return($n);
  }

#------------------------------------------------------------
# create a GenBanl file from flat sequence
#

sub flat2gb ( $ $ $ )
  {
    my($file, $locus, $seq) = @_;
    my($i, $j, $l);

    open(FF,">$file") || die "can't create file $file";
    $l = length($seq);
    print FF "LOCUS\t$locus\n";
    print FF "ORIGIN\n";
    for ($i = 0; $i < $l; $i+=60)
      {
	printf FF "%6d ", $i+1;
	for ($j=0; $j < 6 && $i + $j*10 < $l; $j++)
	  {
	    printf FF " %s", substr($seq,$i + $j*10, 10);
	  }
	print FF "\n";
      }
    print FF "//\n";
    close(FF);
  }


#------------------------------------------------------------
# Reads in an entire (multiple) fasta file and returns a hash in which
# the keys are the identifiers of the sequences (without the '>')
# and the values are the sequences themselves. Created by cesim 14/01/2002
#

sub fasta2hash ( $ )
 {
  my ($file,$key,$value);
  my (%fasta_hash);
  $file=$_[0];
  open (IN,$file);
  while (<IN>)
   {
    chomp;
    if (/^>(\S+)/)
     {
      $key=$1;
     } #if (/^>(\w)$/)
    else 
     {
      $key || die "File $file is not a fasta file!";
      s/\s+//g;
      $fasta_hash{$key}.=$_;
     } #else 
   } #while (<IN>)
  close IN;
  return (%fasta_hash); 
 } #fasta2hash ( $ )

1;
