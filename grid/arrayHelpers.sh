#!/bin/bash

tablerow2array() {
	file=$1
	row=$2
	delimiter=$3
	delimiter=${delimiter:-"\t"}
	line=$(perl -e 'open(F,$ARGV[0]); my $i=0; while (<F>) { my @a=split /$ARGV[2]/; print join(";", @a),"\n" if $i==$ARGV[1]; $i++;}' $file $row $delimiter)
	IFS=";" read -ra array <<< "$line"
	export array
}
containsElement () {
  local e
  for e in "${@:2}"; do [[ "$e" == "$1" ]] && return 0; done
  return 1
}
