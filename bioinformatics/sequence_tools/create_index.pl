#!/usr/bin/env perl
#Generate index files for fasta|embl|swissprot retrieval
#DL 31.03.2004
#USAGE: create_index.pl [fasta|embl|swiss|genbank] targetDBM source(s)

use strict;
use Bio::Index::Fasta;
use Bio::Index::EMBL;
use Bio::Index::Swissprot;
use Bio::Index::GenBank;

use strict;

my $flag = shift;
my $Index_File_Name = shift || die "USAGE $0 [fasta|embl|swiss|genbank] targetDBM source(s)";
my $inx;

if ($flag =~ /fas/i) {
    $inx = Bio::Index::Fasta->new(
				     '-filename' => $Index_File_Name, '-write_flag' => 1);
}
elsif ($flag =~/embl/i) {
    $inx = Bio::Index::EMBL->new('-filename' => $Index_File_Name,
				    '-write_flag' => 'WRITE');
}
elsif ($flag =~ /swiss/i) {
    $inx = Bio::Index::Swissprot->new('-filename' => $Index_File_Name,
					 '-write_flag' => 'WRITE');
}
elsif ($flag =~ /genbank/i) {
    $inx = Bio::Index::GenBank->new('-filename' => $Index_File_Name,
					 '-write_flag' => 'WRITE');
}

$inx->make_index(@ARGV);
