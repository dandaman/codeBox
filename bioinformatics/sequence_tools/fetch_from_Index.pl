#!/usr/bin/env perl
#DL 19.08.2004


use warnings;
use strict;
use Bio::Index::Fasta;
use Bio::Index::GenBank;
use Bio::Index::Swissprot;
use Bio::Index::EMBL;
use Bio::SeqIO;
use Bio::Seq;
use IO::String;
use Getopt::Std;

my $USAGE = "USAGE: $0 -h | -f format-i indexfile ACCESSION NUMBER\n\n";

my %opts;

getopts('i:f:h', \%opts);

die $USAGE if ($opts{'h'});

die $USAGE  unless ($opts{'h'} || ($opts{'i'}));
$opts{f}='fasta' unless $opts{f};
my $index;
$index = Bio::Index::Fasta->new('-filename' => $opts{'i'}) if ($opts{f} eq 'fasta');
$index = Bio::Index::GenBank->new('-filename' => $opts{'i'}) if ($opts{f} eq 'genbank');
$index = Bio::Index::Swissprot->new('-filename' => $opts{'i'}) if ($opts{f} eq 'swiss');
$index = Bio::Index::EMBL->new('-filename' => $opts{'i'}) if ($opts{f} eq 'embl');

my $out =  Bio::SeqIO->new('-format' => $opts{f});
foreach my $acc (@ARGV) {
	#$acc =~ s/\s*|\?//g;
	#$acc = uc($acc);

	die $USAGE unless ($acc);


	my $seq= $index->fetch($acc);

	die "Sequence $acc not found in $opts{'d'}:$opts{'n'}!" unless ($seq);

	$out->write_seq($seq);
}
exit;
