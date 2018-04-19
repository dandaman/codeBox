#!/bin/perl 
use warnings;
use strict;
#use lib qq($ENV{HOME}/bioperl/bioperl-live);
#use lib qq($ENV{HOME}/bioperl/bioperl-run/lib);
use Bio::SeqIO;
use Bio::Factory::EMBOSS;
use Bio::AlignIO;
use File::Temp qw(tempfile);
use Getopt::Long;
use Bio::Align::Utilities qw(aa_to_dna_aln);

my ($help, @model, $USAGE, $KaKs_Calculator, $gapopen, $gapextend);

$USAGE                  = "$0 -h | [-m model] CDS_file\n";
$KaKs_Calculator        = "$ENV{HOME}/software/KaKs_Calculator2.0/src/KaKs_Calculator -i %s -o %s %s"; #inaxt out -m model
$gapopen                = '10.0';
$gapextend              = '0.5';

GetOptions ("help"      => \$help,
            "model=s"   => \@model,
            "open=f"    => \$gapopen,
            "extend=f"  => \$gapextend,
);

@model                  = ("YN") unless @model;


die $USAGE if ($help || ! @ARGV) ;
my $in                  = Bio::SeqIO->new(-format=>"fasta",-file=>$ARGV[0]);
my $factory             = Bio::Factory::EMBOSS ->new()->program("needle");

my (%nt,@aa,%names, $N, %stops);

while (my $s = $in->next_seq) {
    $s->description();
    $N++;
    $names{$N}=$s->display_id;
    $s->display_id($N);
    my $p = $s->translate;
    my $pp= $p->seq;
    $pp=~  s/\*$//;
    $stops{$N}=($pp=~ tr/\*/\*/);
    $nt{$s->display_id}=$s;
    push @aa, $p;
}

my %seen;
for (my $i=0; $i< @aa-1; $i++) {
    for (my $e=0; $e< @aa; $e++) {
		next if $i == $e;
        next if exists $seen{$i}{$e} || exists $seen{$e}{$i};
        my ($fh,$filename)      = tempfile( UNLINK=>1, SUFFIX => '.needle');
        $factory->run({ 
                -asequence      => $aa[$i],
                -bsequence      => $aa[$e],
                -gapopen        => $gapopen,
                -gapextend      => $gapextend,
                -outfile        => $filename});
        my $alnin       = Bio::AlignIO->new(-format=>"emboss", -file=>$filename);
        my $aln         = $alnin->next_aln;
		close($fh);
		undef($fh);
		undef($filename);
	unless ($aln && $aln->isa("Bio::SimpleAlign") ) {
		printf STDERR "No alignment for %s and %s\n", map ({ $names{$_} } ($aa[$i]->display_id, $aa[$e]->display_id));
		$seen{$i}{$e}++;
		$seen{$e}{$i}++;
		next;
	}
        my $dna_aln     = &aa_to_dna_aln($aln, \%nt);
	my $clipped	= $dna_aln->remove_gaps;
	if ($clipped->length>0) { 
		my $axt         = write_axt($clipped);
		my @col         = run_KaKsCalc($axt,$e,\%names, $dna_aln,\%nt, \%stops);
	}
	else {
		printf STDERR "No gap-free alignment left for %s and %s\n" , map { $names{$_->id} } $aln->each_seq;
		my $AA=Bio::AlignIO->new(-format=> "clustalw",-fh=>\*STDERR);
		$AA->write_aln($aln);
	}
        $seen{$i}{$e}++;
        $seen{$e}{$i}++;
    }
}

##SUBS
sub run_KaKsCalc {
    my ($axt,$count,$names, $aln, $nt,$stops) = @_;
    my ($fh,$filename)  = tempfile(UNLINK=>1,SUFFIX=>".KaKs");
    my $cmd             = sprintf($KaKs_Calculator, $axt, $filename, join(" ", map { "-m $_" } @model));
    my $out             = `$cmd`;
    die "Kaks_Calculator failed: $out\n" unless $out =~ /Mission accomplished/;
    my ($c);
    while (<$fh>) {
        chomp;
        my @a = split/\t/;
        if (! $c && $a[0] =~ /^Sequence/) {
            $a[0]= "seq2";
            print join("\t", "seq1",@a, "pid","num_residues","seq1_length","seq2_length", "seq1_stops", "seq2_stops"),"\n" if $count ==1;
        }
        else {
            my @b = map{ $names->{$_} } split/\&/, shift @a;
            print join("\t",@b, @a, $aln->average_percentage_identity, $aln->num_residues, $nt->{1}->length, $nt->{2}->length, $stops->{1}, $stops->{2}),"\n";
        }
        $c++;
    }
	close($fh);
	undef($fh);
	undef($filename);
}

sub write_axt {
    my $aln = shift;
    my ($fh,$filename)      = tempfile(UNLINK=>1,SUFFIX=>".axt");
    my @s = $aln->each_seq;
    print $fh join("&", map { $_->id} @s),"\n";
    print $fh join("\n", map { $_->seq } @s),"\n";
	close($fh);
    return $filename;
}
