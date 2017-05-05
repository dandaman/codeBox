#!/usr/bin/perl -w
#Created by Cedric Simillion on Thu Oct  3 15:35:56 MEST 2002
#Mutilated by Jan Wuyts on 20030331

=head1 Usage
 
> stripalign3.pl < my_alignment.tfa [-vvc] [-g<max_gap_portion>] [-p<percentile>] 
	[-f<pos_save_file>] [-b<blocksize>]

	-g<max_gap_portion> : <max_gap_portion> is the maximum fraction of gaps tolerated per
						  alignment position. (default = 0.2)
	-p<percentile>      : which <percentile> to calculate. (default = 50)
	-f<pos_save_file>   : save translation "new positions" -> "old positions" in file
						  <pos_save_file> (default = pos_save)
	-c					: per codon
	-v					: verbose (use double for more effect)
	-b<blocksize>		: blocksize of output fasta file. (default is all on 1 line)
	-b					: blocksize = 50

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
#####################################################################
# CHANGELOG
#####################################################################
# 20030331: indentation redone (jawuy)
# 20030331: added option -f<pos_save_file> (jawuy)
# 20030401: order of output = order of input (jawuy)
# 20030401: added option -b<blocksize> (jawuy)
# 20030401: added options -v, -vv (jawuy)
#####################################################################
use strict;
use lib qw(/home/lang/scripts/TreePipe/msa);
use bioutils;
use Statistics::Descriptive;

#===============================================================================
# Initialisation
#===============================================================================
my ($per_codon,$c,$savefile,$i,$last_anchor,$stripping,$max_gap_portion,$percentile,$cut,$totpos);
my (@align,@ids,@sequences,$sep,$be_verbose,$blocksize);
my (%align,%blosum,%savep);
$blocksize=0;
$be_verbose=0;
$max_gap_portion = 0.2;
$percentile = 50;
$savefile = "pos_save";
$per_codon = 0;
ARG: while (@ARGV && $ARGV[0] =~ s/^-(?=.)//) {
	OPT: for (shift @ARGV) {
		m/^$/             && do {                        next ARG; };
		m/^-$/            && do {                        last ARG; };
		s/^c//            && do { $per_codon++;          redo OPT; };
		s/^v//            && do { $be_verbose++;         redo OPT; };
		s/^g([0-9.e-]+)// && do { $max_gap_portion = $1; next ARG; };
		s/^p([0-9.e-]+)// && do { $percentile = $1;      next ARG; };
		s/^f(.+)//        && do { $savefile = $1;        next ARG; };
		s/^b([0-9]*)//    && do { $blocksize = $1 || 50; redo OPT; };
		print STDERR "Unknown or incomplete option: $_\n";
	}
}
($per_codon) ? ($per_codon=3) : ($per_codon=1);
$be_verbose-1>0 && print STDERR "per_codon : $per_codon\n";
$be_verbose-1>0 && print STDERR "max_gap_portion : $max_gap_portion\n";
$be_verbose-1>0 && print STDERR "percentile : $percentile\n";
$be_verbose-1>0 && print STDERR "savefile : $savefile\n";
$be_verbose-1>0 && print STDERR "blocksize : $blocksize\n";

%align=&input2hash;
%blosum=&read_blosum;
$be_verbose && print STDERR "Stripping...\n";

#===============================================================================
# Main
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Construct a data structure for editing the alignment
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$c=0;
foreach (keys %align) {
	my ($n);
	$ids[$c]=$_;
	$n=0;
	for ($i=0;$i<length($align{$_});$i+=$per_codon) {
		$align[$n][$c]=substr($align{$_},$i,$per_codon);
		$n++;
	}
	$c++; 
}
#print STDERR ".\n";

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Strip the alignment
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
$last_anchor=0;
$cut=0;
$sep=".";
for ($i=0;$i<scalar(@align);$i++) { # for every site in alignment
	$be_verbose-1>0 && print STDERR $sep.($i+$cut-1);$sep=".";
	if (&gap(@{$align[$i]})|| ($i==0)) {
		$last_anchor=$i;
		while (!&anchor(@{$align[$last_anchor]}) && ($last_anchor>-1)) {$last_anchor--}
		$last_anchor++;
		do {
			splice(@align,$last_anchor,1);
			$cut++;
		} until ( ($last_anchor == scalar(@align)) || &anchor(@{$align[$last_anchor]}) );
		$i=$last_anchor;
		$sep=">";
	}
	$savep{$i}=$i+$cut;
}
	$be_verbose-1>0 && print STDERR "\n";
splice (@align,$last_anchor+1);

$totpos=$i-1;
$be_verbose && print STDERR "Cut out $cut alignment sites\n";
open (OUT,">".$savefile) or die "can't open file $savefile : $!";

for($i=0;$i<=$totpos;$i++) {
	print OUT $i."\t".$savep{$i}."\n";
}
close OUT;
 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Store the edited alignment again in the original hash structure
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%align=();
foreach (@align) {
	for ($i=0;$i<scalar(@{$_});$i++) {
		$align{$ids[$i]}.=$$_[$i];
	}
}
#print STDERR ".\n";

foreach (@sequences) {
	$blocksize && $align{$_} =~ s/(.{$blocksize})/$1\n/g;
	print STDOUT ">$_\n$align{$_}\n";
}

#===============================================================================
# Subroutines
#===============================================================================
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub fasta2hash
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub input2hash {
	my ($key);
	my (%fasta_hash);
	while (<STDIN>) {
		chomp;
		if (/^>(\S+)/) {
			$key=$1;
			push @sequences, $key;
		} else {
			$key || die "Input is not in fasta format!";
			s/\s+//g;
			s/[^ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv-]/X/g;
			$fasta_hash{$key}.=$_;
		}
	} #while (<STDIN>)
	return (%fasta_hash);
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub gap
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub gap {
	my ($test,$gaps);
	$gaps=0;
	$test= '-' x $per_codon;
	foreach (@_) {
		($_ eq $test) && ($gaps++);
	}
	($gaps/scalar(@_)>$max_gap_portion) ? (return -1) : (return 0); 
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub anchor
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
sub anchor {
	my (@row)=@_;
	my ($return,$test,$gap);
	$return=0;
	$gap='-' x $per_codon;
	unless (&gap(@row)) {
		my ($i,$blosums);
		$blosums=Statistics::Descriptive::Full->new();
		if ($per_codon==3) {
			@row=map  {$_=&translate($_)}(@row);
		}
		$test=shift(@row);
		foreach $i (0..$#row-1) {
			my ($j);
			($row[$i] eq $gap) && (next);
			foreach $j ($i+1..$#row) {
				($row[$j] eq $gap) && (next);
				$blosums->add_data($blosum{$row[$i]}{$row[$j]});
			}
		}
		($blosums->percentile($percentile)>=0) && ($return=-1) ;
#		($blosums->median>=0) && ($return=-1);
	}
	return ($return); 
}

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# sub read_blosum
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sub read_blosum {
	my $blosum_file="/home/lang/scripts/TreePipe/msa/blosum62.txt";
	my (@residus);
	my (%matrix);
	open (IN,$blosum_file) || die "Could not open $blosum_file: $!\n";
	(undef,@residus)=split(/[\t\n]/,<IN>);
	while (<IN>) {
		chomp;
		my ($this,@row)=split(/\t/);
		@{$matrix{$this}}{@residus}=@row;
	}
	return (%matrix) 
}
