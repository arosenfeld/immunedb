#!/bin/bash

if [ "$#" -ne 4 ]
then
    echo "Usage: $0 species locus in_dir out_dir"
    exit
fi

species=$1
locus=$2
in_dir=$3
out_dir=$4

mkdir -p $out_dir

for fn in $in_dir/*.fast[aq]
do
    if [[ $locus == T* ]]
    then
        seq_type=TCR
    else
        seq_type=Ig
    fi

    out_fn=$out_dir/`basename $fn .fastq`.airr.tsv
    echo "Running IgBLAST on $fn (saving to $out_fn)"

    cat $fn | ([[ $fn == *.fastq ]] && sed -n '1~4s/^@/>/p;2~4p' || cat -) | \
        ${IGDATA}/bin/igblastn \
        -germline_db_V ${IGDATA}/database/${species}/${locus}V.edited.fasta \
        -germline_db_D ${IGDATA}/database/${species}/${locus}D.edited.fasta \
        -germline_db_J ${IGDATA}/database/${species}/${locus}J.edited.fasta \
        -outfmt 19 \
        -num_threads 12 \
        -domain_system imgt \
        -ig_seqtype $seq_type \
        -auxiliary_data ${IGDATA}/optional_file/${species}_gl.aux \
        -organism ${species} > $out_fn
done
