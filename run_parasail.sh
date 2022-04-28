#!/bin/bash

BASENAME=$1
a=$2
LENGTH=$3
PROGRAM=$4

cat <(echo "@SQ	SN:1	LN:$LENGTH") \
    <(echo "@SQ	SN:2	LN:$LENGTH") \
    <(
paste <(cat $BASENAME.sim.$a.refs.fa | paste - -) \
      <(cat $BASENAME.sim.$a.$PROGRAM.fa | paste - -) \
      <(cat $BASENAME.sim.$a.$PROGRAM.fa | tr -d ">" | paste - - | tr "ATGC" "TACG" | rev | awk -F"\t" '{print ">2\t"$1}') | tr '$' 'N' |
while read F; do
    echo "$F" | tr "\t" "\n" | ../parasail_aligner -a parasail_nw_trace_striped_avx2_256_16 -c 1 -e 1 -o 1 -M 1 -X 1 -x -s 0 -l 0 -i 0 -d -b 64 -g /dev/stdout -O SAM | head -n 2
done) | samtools view -h
