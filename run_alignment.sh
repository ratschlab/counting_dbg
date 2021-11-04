#!/bin/bash

METAGRAPH=../../metagraph

if [[ -z $1 ]]; then
echo "Please provide an input FASTA/Q file"
exit
fi

INFILE=$1
BASENAME=$(basename $INFILE .fa)

mkdir -p swap

echo "Generating Metagraph index"
$METAGRAPH build -v -k 31 -o $BASENAME $INFILE
$METAGRAPH annotate -v -i $BASENAME.dbg --coordinates --anno-filename -o $BASENAME $INFILE
cp $BASENAME.column.annodbg $BASENAME.base.column.annodbg
cp $BASENAME.column.annodbg.coords $BASENAME.base.column.annodbg.coords

$METAGRAPH transform_anno -v --anno-type column_coord --coordinates -o $BASENAME $BASENAME.column.annodbg

$METAGRAPH transform_anno -v --anno-type row_diff --coordinates --row-diff-stage 0 -i $BASENAME.dbg -o out --disk-swap swap $BASENAME.column.annodbg
$METAGRAPH transform_anno -v --anno-type row_diff --coordinates --row-diff-stage 1 -i $BASENAME.dbg -o out --disk-swap swap $BASENAME.column.annodbg
$METAGRAPH transform_anno -v --anno-type row_diff --coordinates --row-diff-stage 2 -i $BASENAME.dbg -o out --disk-swap swap $BASENAME.column.annodbg
$METAGRAPH transform_anno -v --anno-type row_diff_coord -i $BASENAME.dbg -o $BASENAME $BASENAME.column.annodbg
$METAGRAPH transform_anno -v --anno-type row_diff_brwt_coord -i $BASENAME.dbg -o $BASENAME $BASENAME.column.annodbg --greedy --fast --subsample 1000000

echo "Generating minimap2 index"
minimap2 -d $BASENAME.illumina.mmi $INFILE -x sr
minimap2 -d $BASENAME.pacbio.mmi $INFILE -x map-pb

echo "Generating vg index"
vg construct -r $BASENAME.fa -p > $BASENAME.vg
mkdir tmp
vg index -x $BASENAME.xg -g $BASENAME.gcsa $BASENAME.vg -b tmp

echo "Generating BLAST index"
makeblastdb -in $BASENAME.fa -parse_seqids -blastdb_version 5 -title $BASENAME -dbtype nucl -out $BASENAME
makembindex -input $BASENAME

echo "Generating Puffaligner"
../pufferfish index -r $BASENAME.fa --keepFixedFasta -k 31 -o $BASENAME.pufferfish/ >$BASENAME.pufferfish.out 2>$BASENAME.pufferfish.log

echo "Running evaluation"
./evaluate_alignment.sh $BASENAME.fa
