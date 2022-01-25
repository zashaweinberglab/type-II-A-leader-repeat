#!/bin/bash -x

addFile=$1
foldingFlags=$2 # you should put them in quotes e.g., "--extend-repeat-by 15 --no-within-repeat --shuffle-repeat"

datasets=$3 # "II-A--training II-C--training" # take out of dirs -- there's only two, and then it's easier to name the output files
datasetsDesc=$4

outFileList=

for d in $datasets; do
    echo $d
    outFile=out-$d-$addFile.tab;
    python3 /u/homes/zasha/bliss/code/motifs_2007/rpl/ChaseFolding.py $d.tab -o $outFile --num-sample 1000 --num-boltzmann 1000 --print-actual-value --cpu 12 $foldingFlags;
    outFileList="$outFileList $outFile"
done

summarizeDir=summarize-$addFile-$datasetsDesc
mkdir -p $summarizeDir
python3 ~/bliss/code/motifs_2007/rpl/ChaseSummarize.py --out-dir $summarizeDir $outFileList
