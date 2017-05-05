#!/bin/bash
# AUTHOR:   lang
# FILE:     /home/lang/projects/reviews/mads/several/do_codeml.sh
# CREATED:  08:49:37 26/08/2009
# MODIFIED: 08:49:37 26/08/2009
INPUT_NA1=$1

INPUT_NA=$(perl -MFile::Temp -e 'my (@a)= File::Temp::tempfile(-CLEANUP=>1) ; print "$a[1]\n";')
KAKSOUT=$(perl -MFile::Temp -e 'my (@a)= File::Temp::tempfile(-CLEANUP=>1) ; print "$a[1]\n";')
DTVOUT=$(perl -MFile::Temp -e 'my (@a)= File::Temp::tempfile(-CLEANUP=>1) ; print "$a[1]\n";')

perl -p -e 's/>(\S+).*$/>$1/g;' $INPUT_NA1  >  $INPUT_NA

perl pairwiseKaKs.pl $INPUT_NA -m GMYN > $KAKSOUT

perl -e 'while (<>) {if (/^>(\S+)/) {$l{$1}++;}} foreach my $a (sort keys %l) {foreach my $v (sort keys %l) {print "$a\t$v\n" if !($a eq $v) and !($s{$v.$a} or $s{$a.$v}); $s{$a.$v}++; $s{$v.$a}++; }}' $INPUT_NA > $INPUT_NA.pairs

perl $HOME/scripts/TreePipe/msa/4dtv.pl -fasta $INPUT_NA -pairs $INPUT_NA.pairs > $DTVOUT

perl -e 'open(F, $ARGV[0]); while (<F>) { chomp; my @a=split/\t/; my @b=splice(@a,0,2); @a=map {s/^\-+$/NA/; $_} @a; $l{$b[0]}{$b[1]}=\@a;} open(F,$ARGV[1]); while (<F>) { chomp; @a=split/\t/; print join("\t", @a, exists $l{$a[0]}{$a[1]} ?  @{$l{$a[0]}{$a[1]}} : ("NA","NA", "NA")),"\n"; }' $DTVOUT $KAKSOUT

rm $INPUT_NA.pairs $INPUT_NA $INPUT_PROT $KAKSOUT $DTVOUT
