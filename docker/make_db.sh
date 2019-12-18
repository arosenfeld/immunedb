#!/bin/bash

if [ "$#" -lt 2 ]
then
    echo "Usage: $0 out_dir file [file [file ...]]"
    exit
fi

outdir=$1
files=${@:2}

mkdir -p $outdir
for fn in $files
do
    outname=$outdir/`basename $fn .fasta`.gapped.fasta
    sed -E 's/>([^|]+)\|([^|]+)\|.*/>\2/; s/\./-/g' $fn > $outname
    outname=$outdir/`basename $fn .fasta`.edited.fasta
    edit_imgt_file.pl $fn > $outname
    makeblastdb -parse_seqids -dbtype nucl -in $outname
done
