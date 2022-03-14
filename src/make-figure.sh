#!/bin/bash 
set -x
set -e

# example
# ./src/make-figure.sh data/figure-data/II-A-systems-for-figure.tab II-A-figure

# make-figure.sh <input-systems> <outbase-base-name> <flags-for-Python-script>

rm -f temp.* outfile outtree

cat $1 | cut -f 1,3 | sed "s/^/>/" | tr "\t" "\n" > temp.fasta

clustalo --force --in temp.fasta -t RNA --guidetree-out=temp.guidetree --out=temp.sto --outfmt=st --output-order=tree-order --distmat-out=temp.distmat --full

perl src/abbrev-distmat.pl temp.distmat temp.phylip-distmat temp.phylip-tree temp.my-names-phylip-tree

src/make-neighbor-input.sh temp.phylip-distmat | neighbor
rm outfile
mv outtree temp.phylip-tree

perl src/abbrev-distmat.pl temp.distmat temp.phylip-distmat temp.phylip-tree temp.my-names-phylip-tree

python3 src/ChaseFigureFromFolding.py $3 $1 -o $2 --use-prob --no-within-repeat --phylip-tree temp.my-names-phylip-tree

rm -f $2
rm -f $2.html
