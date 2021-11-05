#!/bin/bash

METAGRAPH=../../metagraph

if [[ -z $1 ]]; then
echo "Please provide an input FASTA/Q file"
exit
fi

INFILE=$1
BASENAME=$(basename $INFILE .fa)
FLAGS="-v --align-min-seed-length 19"
MINIMAP2_FLAGS="-a --eqx -t 1"
COORD_TYPE="row_diff"
SEED="1626366265"
ILLUMINA_LENGTH=150
PACBIO_LENGTH=10000
ILLUMINA_XDROP=100
PACBIO_XDROP=27
ILLUMINA_SEED_LENGTH=10000
PACBIO_SEED_LENGTH=19

echo "Simulating Illumina reads"
art_illumina -ss HS25 -i $INFILE -p -c 1280 -l $ILLUMINA_LENGTH -m 200 -s 10 -o $BASENAME.sim -rs $SEED -sam
reformat.sh in1=$BASENAME.sim1.fq in2=$BASENAME.sim2.fq out=$BASENAME.sim.illumina.fq overwrite=true fastawrap=-1
cat <(samtools view -H $BASENAME.sim.sam | awk '$1=="@SQ"{print $1"\t"$2"\t"$(NF); next}{print}') <(grep -v "^@" $BASENAME.sim.sam) | samtools view -h -b > $BASENAME.sim.bam
bedtools getfasta -fi $BASENAME.fa -bed <(bedtools bamtobed -i $BASENAME.sim.bam) > $BASENAME.sim.illumina.refs.fa
reformat.sh in=$BASENAME.sim.illumina.fq out=$BASENAME.sim.illumina.fa fastawrap=-1 overwrite=true

echo "Simulating PacBio reads"
PBSIM_DEPTH=0.0504 # chr22
#PBSIM_DEPTH=0.4284 # E. coli
../pbsim --data-type CLR --depth $PBSIM_DEPTH --model_qc ../model_qc_clr --data-type CLR --length-mean $PACBIO_LENGTH --length-sd 0 --length-min $PACBIO_LENGTH --length-max $PACBIO_LENGTH $BASENAME.fa --prefix $BASENAME --seed $SEED
cat ${BASENAME}_0001.fastq | paste - - - - | awk '!($2~/N/)' | shuf --random-source=../seed | tr "\t" "\n" > $BASENAME.sim.pacbio.fq
cat ${BASENAME}_0001.maf | paste - - - - | awk -F"\t" '!($(NF-1)~/N/)' | shuf --random-source=../seed | tr "\t" "\n" > $BASENAME.sim.maf
cat $BASENAME.sim.maf | paste - - - - | awk '{print ">1\n"$(NF-7)}' | tr -d "-" > $BASENAME.sim.pacbio.refs.fa
reformat.sh in=$BASENAME.sim.pacbio.fq out=$BASENAME.sim.pacbio.fa fastawrap=-1 overwrite=true

for a in "illumina" "pacbio"; do
if [[ $a == "illumina" ]]; then
    LENGTH=$ILLUMINA_LENGTH
    XDROP=$ILLUMINA_XDROP
    SEED_LENGTH=$ILLUMINA_SEED_LENGTH
else
    LENGTH=$PACBIO_LENGTH
    XDROP=$PACBIO_XDROP
    SEED_LENGTH=$PACBIO_SEED_LENGTH
fi

echo "Aligning to whole graph"
/usr/bin/time -v $METAGRAPH align $FLAGS --align-xdrop $XDROP -i $BASENAME.dbg $BASENAME.sim.$a.fq > $BASENAME.sim.$a.metagraph_base.out 2> $BASENAME.sim.$a.metagraph_base.log
awk -F"\t" '{print ">1\n"$4}' $BASENAME.sim.$a.metagraph_base.out > $BASENAME.sim.$a.metagraph_base.out.fa
echo "   parasail"
../run_parasail.sh $BASENAME $a $LENGTH metagraph_base.out > $BASENAME.sim.$a.metagraph_base.out.results.sam

echo "Aligning by coords"
/usr/bin/time -v $METAGRAPH align $FLAGS --align-chain --align-max-seed-length $SEED_LENGTH --align-xdrop 100 -i $BASENAME.dbg $BASENAME.sim.$a.fq -a $BASENAME.${COORD_TYPE}_coord.annodbg > $BASENAME.sim.$a.metagraph_coord.out 2> $BASENAME.sim.$a.metagraph_coord.log
awk -F"\t" '{print ">1\n"$4}' $BASENAME.sim.$a.metagraph_coord.out > $BASENAME.sim.$a.metagraph_coord.out.fa
echo "   parasail"
../run_parasail.sh $BASENAME $a $LENGTH metagraph_coord.out > $BASENAME.sim.$a.metagraph_coord.out.results.sam

echo "Aligning with minimap2"
/usr/bin/time -v minimap2 $MINIMAP2_FLAGS $BASENAME.$a.mmi $BASENAME.sim.$a.fq > $BASENAME.sim.$a.minimap2.sam 2> $BASENAME.sim.$a.minimap2.log
samtools view -h -F 256 $BASENAME.sim.$a.minimap2.sam | awk -F"\t" -v OFS='\t' '$3=="*"{$2="0"; $3="'$BASENAME'"; $4="1"; $6="'$LENGTH'X"}{print}' > $BASENAME.sim.$a.minimap2.nosec.sam
bedtools getfasta -fi $BASENAME.fa -bed <(bedtools bamtobed -i $BASENAME.sim.$a.minimap2.nosec.sam) | paste - - | awk '{print ">1\n"$2}' > $BASENAME.sim.$a.minimap2.sam.fa
echo "   parasail"
../run_parasail.sh $BASENAME $a $LENGTH minimap2.sam > $BASENAME.sim.$a.minimap2.sam.results.sam

echo "Aligning with vg"
/usr/bin/time -v vg map -d $BASENAME -f $BASENAME.sim.$a.fq -t 1 --surject-to bam 2> $BASENAME.sim.$a.vg.log | samtools view -h > $BASENAME.sim.$a.vg.sam
bedtools getfasta -fi $BASENAME.fa -bed <(bedtools bamtobed -i $BASENAME.sim.$a.vg.sam) | paste - - | awk '{print ">1\n"$2}' > $BASENAME.sim.$a.vg.sam.fa
echo "   parasail"
../run_parasail.sh $BASENAME $a $LENGTH vg.sam > $BASENAME.sim.$a.vg.sam.results.sam

echo "Aligning with BLAST"
/usr/bin/time -v blastn -query $BASENAME.sim.$a.fa -db $BASENAME -use_index true -max_hsps 1 -outfmt "6 qseqid sseq" 2> $BASENAME.sim.$a.blast.log | tr -d "-" > $BASENAME.sim.$a.blast.out
python3 ../parse_blast.py $BASENAME.sim.$a.blast.out $BASENAME.sim.$a.fa $BASENAME.sim.$a.refs.fa > $BASENAME.sim.$a.blast.out.fa
echo "   parasail"
../run_parasail.sh $BASENAME $a $LENGTH blast.out > $BASENAME.sim.$a.blast.out.results.sam

echo "Aligning with Pufferfish"
/usr/bin/time -v ../pufferfish align -i $BASENAME.pufferfish/ --read $BASENAME.sim.$a.fq --genomicReads --primaryAlignment -o $BASENAME.sim.$a.pufferfish.sam > $BASENAME.sim.$a.pufferfish.log 2>&1
bedtools getfasta -fi $BASENAME.fa -bed <(bedtools bamtobed -i $BASENAME.sim.$a.pufferfish.sam) | paste - - | awk '{print ">1\n"$2}' > $BASENAME.sim.$a.pufferfish.sam.fa
echo "   parasail"
../run_parasail.sh $BASENAME $a $LENGTH pufferfish.sam > $BASENAME.sim.$a.pufferfish.sam.results.sam

done

