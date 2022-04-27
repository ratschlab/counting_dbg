## Prepare and build graph

```bash
KMC=~/projects/projects2014-metagenome/metagraph/build_release/KMC/kmc;
DIR=~/metagenome/data/kingsford_21;
# rm -r $DIR/kmc_21_filtered_2gb;
mkdir $DIR/kmc_21_filtered_2gb;
mkdir $DIR/logs;

for cutoff in {1,3,10,20,50}; do
  ids=$DIR/kingsford_${cutoff}.txt;
  bsub -J "filter[1-$(cat $ids | wc -l)]%800" \
       -o $DIR/logs/kmc_count_2gb.lsf \
       -W 4:00 \
       -n 1 -R "rusage[mem=5000] span[hosts=1] select[model==XeonGold_6140]" \
            "id=\\\$(sed -n \${LSB_JOBINDEX}p $ids); \
            mkdir ~/metagenome/scratch/nobackup/stripe_1/\\\${id}.kmc_cache; \
            file=~/metagenome/raw_data/kingsford/data_fasta/\\\${id}.fasta.gz; \
            /usr/bin/time -v $KMC -k21 -m1 -sm -ci$cutoff -fm -t2 \
                \\\${file} \
                $DIR/kmc_21_filtered_2gb/\\\$id \
                ~/metagenome/scratch/nobackup/stripe_1/\\\${id}.kmc_cache \
            2>&1 | tee $DIR/kmc_21_filtered_2gb/\\\${id}.log;
            rm -r ~/metagenome/scratch/nobackup/stripe_1/\\\${id}.kmc_cache"
done

mkdir $DIR/unitigs;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "build_single[1-2652]%500" \
     -w "filter*" \
     -o $DIR/logs/build_single.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=20000] span[hosts=1] select[model==XeonGold_6140]" \
        "id=\\\$(sed -n \${LSB_JOBINDEX}p $DIR/kingsford.txt); \
        file=$DIR/kmc_21_filtered/\\\${id}.kmc_suf; \
        /usr/bin/time -v $METAGRAPH build \
            -k 21 \
            --mode canonical \
            --mem-cap-gb 8 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 2 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            \\\$file; \
        /usr/bin/time -v $METAGRAPH transform \
            --to-fasta --primary-kmers \
            -p 2 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg; \
        rm $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg*"

bsub -J "build_graph" \
     -w "build_single" \
     -oo $DIR/logs/build_graph.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=3000] span[hosts=1] select[model==XeonGold_6140]" \
    "find $DIR/unitigs -name \"*.fasta.gz\" \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 21 \
            --mode canonical \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/kingsford_canonical; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --to-fasta --primary-kmers \
            -o $DIR/kingsford_primary \
            $DIR/kingsford_canonical.dbg \
            -p 36; \
    rm $DIR/kingsford_canonical.dbg; \
    /usr/bin/time -v $METAGRAPH build -v \
            -k 21 \
            --mode primary \
            --mem-cap-gb 50 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 72 \
            -o $DIR/kingsford \
            $DIR/kingsford_primary.fasta.gz; \
    rm $DIR/kingsford_primary.fasta.gz; \
    /usr/bin/time -v $METAGRAPH transform -v \
            --state small \
            -o $DIR/kingsford_small \
            $DIR/kingsford.dbg \
            -p 72"
```


## Index with k-mer counts

```bash
WINDOW_SIZE=1;
# WINDOW_SIZE=1000000000;
DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
mkdir $DIR;
mkdir $DIR/logs;
mkdir $DIR/unitigs;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "build_single_${WINDOW_SIZE}[1-2652]%500" \
     -w "filter*" \
     -o $DIR/logs/build_single.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=20000] span[hosts=1] select[model==XeonGold_6140]" \
        "id=\\\$(sed -n \${LSB_JOBINDEX}p $DIR/../kingsford.txt); \
        file=$DIR/../kmc_21_filtered/\\\${id}.kmc_suf; \
        /usr/bin/time -v $METAGRAPH build -v \
            -k 21 \
            --mode canonical \
            --count-kmers --count-width 32 \
            --mem-cap-gb 8 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -p 2 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            \\\$file; \
        /usr/bin/time -v $METAGRAPH clean -v \
            --to-fasta --primary-kmers \
            --smoothing-window ${WINDOW_SIZE} \
            -p 2 \
            -o $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}) \
            $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg; \
        rm $DIR/unitigs/\\\$(basename \\\${file%.kmc_suf}).dbg*"

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
bsub -J "split_${WINDOW_SIZE}" \
     -w "build_single_${WINDOW_SIZE}" \
     -o /dev/null -W 1:00 -n 1 -R "rusage[mem=1000]" \
        "cd $DIR; \
        mkdir -p batches; \
        cd batches; \
        split -d -n r/10 <(find $DIR/unitigs -name "*.fasta.gz" | shuf); \
        mkdir -p ${DIR}/columns;";

DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
for N in {0..9}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "annotate_${WINDOW_SIZE}_${list}" \
         -w "split_${WINDOW_SIZE} && build_graph" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=1500] span[hosts=1] select[model==XeonGold_6140]" \
        "cat $DIR/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i $DIR/../kingsford.dbg \
                --anno-filename \
                --separately \
                --count-kmers --count-width 32 \
                -o ${DIR}/columns \
                -p 36"; \
done

WINDOW_SIZE=1;
# WINDOW_SIZE=1000000000;

# git checkout 7a9027fa8c6c29742c7885f77f90d414c39c5b53
DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE}_new_enc;
# rm -r $DIR;
mkdir $DIR;
mkdir $DIR/rd;
mkdir $DIR/logs;
mkdir $DIR/rd/rd_columns;
ln -s ~/metagenome/data/kingsford_21/kingsford.dbg ${DIR}/rd/graph.dbg;


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "count_${WINDOW_SIZE}_rd_brwt" \
     -oo ${DIR}/logs/count_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
    "find ${DIR}/../smoothing_${WINDOW_SIZE}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72; \
    find ${DIR}/../smoothing_${WINDOW_SIZE}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72; \
    find ${DIR}/../smoothing_${WINDOW_SIZE}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff --count-kmers \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72; \
    find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy --fast --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10; \
    /usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 72 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_int_brwt.annodbg";
```

### Query

```bash
DIR=~/metagenome/data/kingsford_21;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
bsub -J "kingsford_query_old" \
     -oo ${DIR}/logs/query_rd_brwt_old.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=8000] span[hosts=1] select[model==XeonGold_6140]" \
    "/usr/bin/time -v $METAGRAPH query --count-labels --fast -v \
            --discovery-fraction 0 \
            -i ~/metagenome/finished_projects/counting_dbg/kingsford_21/kingsford_small.dbg \
            -a ~/metagenome/finished_projects/counting_dbg/kingsford_21/annotation_old.relaxed.row_diff_brwt.annodbg \
            ~/projects/projects2014-metagenome/metagraph/tests/data/transcripts_100.fa";

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_query" \
     -oo ${DIR}/logs/query_rd_brwt.lsf \
     -W 4:00 \
     -n 36 -R "rusage[mem=8000] span[hosts=1] select[model==XeonGold_6140]" \
    "/usr/bin/time -v $METAGRAPH query --count-labels --fast -v \
            --discovery-fraction 0 \
            -i ${DIR}/kingsford_small.dbg \
            -a ${DIR}/annotation.relaxed.row_diff_brwt.annodbg \
            ~/projects/projects2014-metagenome/metagraph/tests/data/transcripts_100.fa";

for WINDOW_SIZE in {1,1000000000}; do
    DIR=~/metagenome/data/kingsford_21/smoothing_${WINDOW_SIZE};
    bsub -J "kingsford_count_query" \
         -oo ${DIR}/logs/query_count_rd_brwt.lsf \
         -W 4:00 \
         -n 36 -R "rusage[mem=2000] span[hosts=1] select[model==XeonGold_6140]" \
        "/usr/bin/time -v $METAGRAPH query --count-labels --query-counts --fast -v \
                --discovery-fraction 0 \
                -i ${DIR}/../kingsford_small.dbg \
                -a ${DIR}/annotation.relaxed.row_diff_int_brwt.annodbg \
                ~/projects/projects2014-metagenome/metagraph/tests/data/transcripts_100.fa";
done
```


## Binary annotation without k-mer counts

```bash
DIR=~/metagenome/data/kingsford_21;

mkdir ${DIR}/columns;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
for N in {0..9}; do
    N=$(printf "%02d" $N);
    list=x$N;
    bsub -J "annotate_${list}" \
         -w "split_1 && build_graph" \
         -oo ${DIR}/logs/annotate_${list}.lsf \
         -W 4:00 \
         -n 18 -R "rusage[mem=1500] span[hosts=1] select[model==XeonGold_6140]" \
        "cat $DIR/smoothing_1/batches/${list} \
            | /usr/bin/time -v $METAGRAPH annotate \
                -i $DIR/kingsford.dbg \
                --anno-filename \
                --separately \
                -o ${DIR}/columns \
                -p 36"; \
done

DIR=~/metagenome/data/kingsford_21;
mkdir $DIR;
mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s ~/metagenome/data/kingsford_21/kingsford.dbg ${DIR}/rd/graph.dbg;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "kingsford_rd_brwt" \
     -oo ${DIR}/logs/rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72; \
    find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72; \
    find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72; \
    find ${DIR}/rd/rd_columns -name \"*.row_diff.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --fast --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10; \
    /usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 72 \
            --relax-arity 32 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_brwt.annodbg";
```


## RowDiff 1.0

```bash
git checkout 0d9feb76a9840b92031c25571cbc0f23ffd1cbe2

DIR=~/metagenome/data/kingsford_21;
mkdir $DIR/rd_old;
mkdir $DIR/rd_old/rd_columns;
ln -s ~/metagenome/data/kingsford_21/kingsford.dbg ${DIR}/rd_old/graph.dbg;

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_1" \
     -w "kingsford_annotate_*" \
     -oo ${DIR}/logs/old_rd_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --max-path-length 100 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd_old/graph.dbg \
            -o ${DIR}/rd_old/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/old_rd_1.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_2" \
     -w "kingsford_rd_old_1" \
     -oo ${DIR}/logs/old_rd_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/smoothing_1/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --optimize \
            --max-path-length 100 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd_old/graph.dbg \
            -o ${DIR}/rd_old/rd_columns/out \
            -p 72 \
            2>&1 | tee ${DIR}/logs/old_rd_2.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_brwt" \
     -w "kingsford_rd_old_2" \
     -oo ${DIR}/logs/old_rd_brwt.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=4000] span[hosts=1]" \
    "find ${DIR}/rd_old/rd_columns -name \"*.row_diff.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt \
            --greedy --fast --subsample 1000000 \
            -i ${DIR}/rd_old/graph.dbg \
            -o ${DIR}/annotation_old \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/old_rd_brwt.log";

DIR=~/metagenome/data/kingsford_21;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test2/metagraph;
bsub -J "kingsford_rd_old_brwt_relax" \
     -w "kingsford_rd_old_brwt" \
     -oo ${DIR}/logs/old_rd_brwt_relax.lsf \
     -W 24:00 \
     -n 12 -R "rusage[mem=5000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            -p 24 \
            --relax-arity 32 \
            -o ${DIR}/annotation_old.relaxed \
            ${DIR}/annotation_old.row_diff_brwt.annodbg \
            2>&1 | tee ${DIR}/logs/old_rd_brwt_relax.log";
```


## HiFi Viruses Index with k-mer coordinates (lossless read compression)

### Compress with tools for seq compression
```bash
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data;
mkdir $DIR/compressed;
mkdir $DIR/compressed/logs;
find ~/metagenome/raw_data/hifi_sra/viruses_hifi_data/ -name "*.fastq.gz" > $DIR/list.txt;
bsub -J "noheader_gzip[1-4132]" \
     -o $DIR/compressed/logs/noheader_gzip.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir $DIR/compressed/\\\${id:0:6}; \
              /usr/bin/time -v zcat \\\$file | sed -n '1~4s/^@/>/p;2~4p' | sed 's/^>.*/>/' | gzip -9 > $DIR/compressed/\\\${id:0:6}/\\\${id}_no_header.fasta.gz ; \
              /usr/bin/time -v spring -l -c -g --fasta-input -t 2 \
                                      -i $DIR/compressed/\\\${id:0:6}/\\\${id}_no_header.fasta.gz \
                                      -w ${DIR}/compressed/\\\${id:0:6} \
                                      -o ${DIR}/compressed/\\\${id:0:6}/\\\${id}_no_header.spring; \
              /usr/bin/time -v zcat \\\$file | paste - - - - | cut -f2 | tr -d '\n' | wc -c > $DIR/compressed/\\\${id:0:6}/\\\${id}_num_bp ; \
          done"


bsub -J "spring[1-4132]" \
     -o $DIR/compressed/logs/spring.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir $DIR/compressed/\\\${id:0:6}; \
              mkdir $DIR/compressed/\\\${id:0:6}/\\\${id}_temp;
              /usr/bin/time -v spring -l -c -g --fasta-input -t 2 \
                                      -i $DIR/compressed/\\\${id:0:6}/\\\${id}_no_header.fasta.gz \
                                      -w ${DIR}/compressed/\\\${id:0:6}/\\\${id}_temp \
                                      -o ${DIR}/compressed/\\\${id:0:6}/\\\${id}_no_header.spring; \
              rm -r $DIR/compressed/\\\${id:0:6}/\\\${id}_temp;
              /usr/bin/time -v zcat \\\$file | paste - - - - | cut -f2 | tr -d '\n' | wc -c > $DIR/compressed/\\\${id:0:6}/\\\${id}_num_bp ; \
          done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data;
mkdir $DIR/blast;
mkdir $DIR/blast/logs;

bsub -J "blast[1-4132]" \
     -o $DIR/blast/logs/build_database.lsf \
     -W 24:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir $DIR/blast/\\\${id:0:6}; \
              mkdir $DIR/blast/\\\${id:0:6}/\\\${id}; \
              cd $DIR/blast/\\\${id:0:6}/\\\${id}; \
              cp $DIR/compressed/\\\${id:0:6}/\\\${id}_no_header.fasta.gz ./; \
              gunzip \\\${id}_no_header.fasta.gz; \
              makeblastdb -in \\\${id}_no_header.fasta -dbtype nucl -out \\\${id}; \
              rm \\\${id}_no_header.fasta; \
          done"


DIR=~/metagenome/data/hifi_sra/viruses_hifi_data;
mkdir $DIR/pufferfish_sparse;
mkdir $DIR/pufferfish_sparse/logs;

bsub -J "pufferfish_sparse[1-4132]" \
     -o $DIR/pufferfish_sparse/logs/build_index.lsf \
     -W 24:00 \
     -n 4 -R "rusage[mem=19000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir $DIR/pufferfish_sparse/\\\${id:0:6}; \
              mkdir $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              cd $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              ~/pufferfish index -r \\\$file -k 31 -s -o \\\${id}; \
          done"
# re-running those going out of RAM
bsub -J "pufferfish_sparse[86,104,3488,3828,3509]" \
     -o $DIR/pufferfish_sparse/logs/build_index_2.lsf \
     -W 48:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir -p $DIR/pufferfish_sparse/\\\${id:0:6}; \
              rm -rf $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              mkdir $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              cd $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              ~/pufferfish index -r \\\$file -k 31 -s -o \\\${id}; \
          done"
# re-running those going out of time
bsub -J "pufferfish_sparse[1354]" \
     -o $DIR/pufferfish_sparse/logs/build_index_3.lsf \
     -W 72:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 14)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir -p $DIR/pufferfish_sparse/\\\${id:0:6}; \
              rm -rf $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              mkdir $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              cd $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              ~/pufferfish index -p 36 -r \\\$file -k 31 -s -o \\\${id}; \
          done"
bsub -J "pufferfish_sparse[1354]" \
     -o $DIR/pufferfish_sparse/logs/build_index_4.lsf \
     -W 24:00 \
     -n 4 -R "rusage[mem=19000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 16)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir -p $DIR/pufferfish_sparse/\\\${id:0:6}; \
              rm -rf $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              mkdir $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              cd $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              ~/pufferfish index -r \\\$file -k 31 -s -o \\\${id}; \
          done"
bsub -J "pufferfish_sparse[1354]" \
     -o $DIR/pufferfish_sparse/logs/build_index_5.lsf \
     -W 24:00 \
     -n 4 -R "rusage[mem=19000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 26)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir -p $DIR/pufferfish_sparse/\\\${id:0:6}; \
              rm -rf $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              mkdir $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              cd $DIR/pufferfish_sparse/\\\${id:0:6}/\\\${id}; \
              ~/pufferfish index -r \\\$file -k 31 -s -o \\\${id}; \
          done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data;
mkdir $DIR/megablast;
mkdir $DIR/megablast/logs;

bsub -J "megablast[1-4132]" \
     -w "blast[*]" \
     -o $DIR/megablast/logs/build_database.lsf \
     -W 24:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
          "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
              file=\\\$(sed -n \\\${i}p $DIR/list.txt); \
              id=\\\$(basename \\\${file%.fastq.gz}); \
              mkdir $DIR/megablast/\\\${id:0:6}; \
              mkdir $DIR/megablast/\\\${id:0:6}/\\\${id}; \
              ln -s $DIR/blast/\\\${id:0:6}/\\\${id}/* $DIR/megablast/\\\${id:0:6}/\\\${id}/; \
              makembindex -input $DIR/megablast/\\\${id:0:6}/\\\${id}/\\\${id}; \
          done"
```

### with Metagraph
```bash
K=31
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data/mtg;
DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/mtg_fork_opt;
mkdir $DIR;
mkdir $DIR/logs;

list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5_test/metagraph;
bsub -J "build_${K}_[1-4132]" \
     -W 24:00 \
     -n 1 -R "rusage[mem=20000] span[hosts=1]" \
    "for i in \\\$(seq \\\$((\\\${LSB_JOBINDEX} * 37 - 37 + 1)) \\\$((\\\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir ${DIR}/\\\${id:0:6}; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}/logs; \
        /usr/bin/time -v $METAGRAPH transform -v \
                --index-ranges 1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph \
                ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg; \
        /usr/bin/time -v $METAGRAPH transform -v \
                --state small \
                --index-ranges 1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph_small \
                ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg; \
    done"


METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5_test/metagraph;
NUM_THREADS=18;
bsub -J "build_${K}_[1-4132]" \
     -o ${DIR}/logs/construct_${K}.lsf \
     -W 24:00 \
     -n 18 -R "rusage[mem=10000] span[hosts=1]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir ${DIR}/\\\${id:0:6}/; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}/logs; \

        /usr/bin/time -v $METAGRAPH build -v -p $NUM_THREADS \
                -k ${K} \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph \
                \\\$file \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/build.log; \

        /usr/bin/time -v $METAGRAPH transform -v -p $NUM_THREADS \
                --state small \
                --index-ranges 1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph_small \
                ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/transform.log; \

        /usr/bin/time -v $METAGRAPH annotate -v -p $NUM_THREADS \
                --coordinates \
                --anno-filename \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/annotation \
                \\\$file \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/annotate.log; \

        mkdir ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 0 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_1.log; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 1 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_2.log; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 2 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_3.log; \

        rm ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg.coords; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.row_count; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.row_reduction; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.succ; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.succ_boundary; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.pred; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.pred_boundary; \

        /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_coord \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/annotation \
                ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.column.annodbg; \
    done"

NUM_THREADS=36;
bsub -J "build_${K}_rerun_[1-4132]" \
     -w "exit(build_${K}_[*])" \
     -o ${DIR}/logs/construct_${K}_rerun.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=10000] span[hosts=1]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir ${DIR}/\\\${id:0:6}/; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}/logs; \

        /usr/bin/time -v $METAGRAPH build -v -p $NUM_THREADS \
                -k ${K} \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph \
                \\\$file \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/build.log; \

        /usr/bin/time -v $METAGRAPH transform -v -p $NUM_THREADS \
                --state small \
                --index-ranges 1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph_small \
                ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/transform.log; \

        /usr/bin/time -v $METAGRAPH annotate -v -p $NUM_THREADS \
                --coordinates \
                --anno-filename \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/annotation \
                \\\$file \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/annotate.log; \

        mkdir ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 0 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_1.log; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 1 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_2.log; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 2 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_3.log; \

        rm ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg.coords; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.row_count; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.row_reduction; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.succ; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.succ_boundary; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.pred; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.pred_boundary; \

        /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_coord \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/annotation \
                ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.column.annodbg; \
    done"


bsub -J "build_${K}_rerun2_[1-4132]" \
     -w "exit(build_${K}_[*]) && exit(build_${K}_rerun_[*])" \
     -o ${DIR}/logs/construct_${K}_rerun2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=40000] span[hosts=1]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir ${DIR}/\\\${id:0:6}/; \
        rm -r ${DIR}/\\\${id:0:6}/\\\${id}; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}; \
        mkdir ${DIR}/\\\${id:0:6}/\\\${id}/logs; \

        /usr/bin/time -v $METAGRAPH build -v -p $NUM_THREADS \
                -k ${K} \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph \
                \\\$file \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/build.log; \

        /usr/bin/time -v $METAGRAPH transform -v -p $NUM_THREADS \
                --state small \
                --index-ranges 1 \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/graph_small \
                ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/transform.log; \

        /usr/bin/time -v $METAGRAPH annotate -v -p $NUM_THREADS \
                --coordinates \
                --anno-filename \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/annotation \
                \\\$file \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/annotate.log; \

        mkdir ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 0 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_1.log; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 1 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_2.log; \

        /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 2 \
                --mem-cap-gb 200 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/out \
                ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg \
                2> ${DIR}/\\\${id:0:6}/\\\${id}/logs/coord_rd_3.log; \

        rm ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/annotation.column.annodbg.coords; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.row_count; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.row_reduction; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.succ; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.succ_boundary; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.pred; \
        rm ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg.pred_boundary; \

        /usr/bin/time -v $METAGRAPH transform_anno -v \
                --anno-type row_diff_coord \
                -i ${DIR}/\\\${id:0:6}/\\\${id}/graph.dbg \
                -o ${DIR}/\\\${id:0:6}/\\\${id}/annotation \
                ${DIR}/\\\${id:0:6}/\\\${id}/rd_columns/annotation.column.annodbg; \
    done"


### Query
K=31
DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/mtg;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;

bsub -J "query_[1-4132]" \
     -o ${DIR}/logs/delta_variants_query.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=40000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \

        mkdir /scratch/mtg_\\\${id}; \
        cp ${DIR}/\\\${id:0:6}/\\\${id}/graph_small.dbg /scratch/mtg_\\\${id}/;
        cp ${DIR}/\\\${id:0:6}/\\\${id}/annotation.row_diff_coord.annodbg /scratch/mtg_\\\${id}/;
        echo \\\${id} \\\$(/usr/bin/time -v $METAGRAPH query -v \
                --query-coords \
                -i /scratch/mtg_\\\${id}/graph_small.dbg \
                -a /scratch/mtg_\\\${id}/annotation.row_diff_coord.annodbg \
                ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                2>&1 1>/dev/null | sed '8q;d' | cut -d' ' -f8) >> ${DIR}/delta_variants_query.times; \
        rm -r /scratch/mtg_\\\${id}; \
    done"

K=31
DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/mtg;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;

bsub -J "align_[1-4132]" \
     -o ${DIR}/logs/delta_variants_align.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=200000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \

        mkdir /scratch/mtg_\\\${id}; \
        cp ${DIR}/\\\${id:0:6}/\\\${id}/graph_small.dbg /scratch/mtg_\\\${id}/;
        cp ${DIR}/\\\${id:0:6}/\\\${id}/annotation.row_diff_coord.annodbg /scratch/mtg_\\\${id}/;
        echo \\\${id} \\\$(/usr/bin/time -v $METAGRAPH align -v \
                --align-chain --align-max-seed-length 19 \
                -i /scratch/mtg_\\\${id}/graph_small.dbg \
                -a /scratch/mtg_\\\${id}/annotation.row_diff_coord.annodbg \
                ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                2>&1 1>/dev/null | tail -n 19 | sed '1q;d' | cut -d' ' -f8) >> ${DIR}/delta_variants_align.times; \
        rm -r /scratch/mtg_\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/megablast;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "query_[1-4132]" \
     -o ${DIR}/logs/delta_variants_query.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/\\\${id}; \
        cp -H ${DIR}/\\\${id:0:6}/\\\${id}/* /scratch/\\\${id}/; \
        echo \\\${id} \\\$(/usr/bin/time -v blastn \
                -query ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -db /scratch/\\\${id}/\\\${id} \
                -use_index true \
                -max_hsps 1 -outfmt '6 qseqid qseq sseq bitscore nident btop' \
                2>&1 1>/dev/null | sed '5q;d' | cut -d' ' -f8) >> ${DIR}/delta_variants_query.times; \
        rm -r /scratch/\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/blast;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "query_[1-4132]" \
     -o ${DIR}/logs/delta_variants_query.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/blast_\\\${id}; \
        cp -H ${DIR}/\\\${id:0:6}/\\\${id}/* /scratch/blast_\\\${id}/; \
        echo \\\${id} \\\$(/usr/bin/time -v blastn \
                -query ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -db /scratch/blast_\\\${id}/\\\${id} \
                -max_hsps 1 -outfmt '6 qseqid qseq sseq bitscore nident btop' \
                2>&1 1>/dev/null | sed '5q;d' | cut -d' ' -f8) >> ${DIR}/delta_variants_query.times; \
        rm -r /scratch/blast_\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/pufferfish_sparse;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "query_[1-4132]" \
     -o ${DIR}/logs/delta_variants_query.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=200000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/pf_\\\${id}; \
        cp -r ${DIR}/\\\${id:0:6}/\\\${id}/\\\${id} /scratch/pf_\\\${id}; \
        echo \\\${id} \\\$(/usr/bin/time -v ~/pufferfish lookup \
                -r ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -i /scratch/pf_\\\${id}/\\\${id}/ \
                2>&1 | tail -n 19 | sed '1q;d' | cut -d' ' -f8) >> ${DIR}/delta_variants_query.times; \
        rm -r /scratch/pf_\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/pufferfish_sparse;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "align_[1-4132]" \
     -o ${DIR}/logs/delta_variants_align.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=200000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/pf_align_\\\${id}; \
        cp -r ${DIR}/\\\${id:0:6}/\\\${id}/\\\${id} /scratch/pf_align_\\\${id}; \
        echo \\\${id} \\\$(/usr/bin/time -v ~/pufferfish align \
                --read ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -i /scratch/pf_align_\\\${id}/\\\${id}/ \
                --genomicReads --primaryAlignment -o /dev/null \
                2>&1 | tail -n 19 | sed '1q;d' | cut -d' ' -f8) >> ${DIR}/delta_variants_align.times; \
        rm -r /scratch/pf_align_\\\${id}; \
    done"


# Test

K=31
DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/mtg;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;

bsub -J "query_[1-1]" \
     -o ${DIR}/logs/delta_variants_query_1.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=40000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \

        mkdir /scratch/mtg_\\\${id}; \
        cp ${DIR}/\\\${id:0:6}/\\\${id}/graph_small.dbg /scratch/mtg_\\\${id}/;
        cp ${DIR}/\\\${id:0:6}/\\\${id}/annotation.row_diff_coord.annodbg /scratch/mtg_\\\${id}/;
        /usr/bin/time -v $METAGRAPH query -v \
                --query-coords \
                -i /scratch/mtg_\\\${id}/graph_small.dbg \
                -a /scratch/mtg_\\\${id}/annotation.row_diff_coord.annodbg \
                ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                2>&1; \
        rm -r /scratch/mtg_\\\${id}; \
    done"


K=31
DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/mtg;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;

bsub -J "align_[1-1]" \
     -o ${DIR}/logs/delta_variants_align_1_chain_optim.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=40000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \

        mkdir /scratch/mtg_\\\${id}; \
        cp ${DIR}/\\\${id:0:6}/\\\${id}/graph_small.dbg /scratch/mtg_\\\${id}/;
        cp ${DIR}/\\\${id:0:6}/\\\${id}/annotation.row_diff_coord.annodbg /scratch/mtg_\\\${id}/;
        /usr/bin/time -v $METAGRAPH align -v \
                --align-chain --align-max-seed-length 19 \
                -i /scratch/mtg_\\\${id}/graph_small.dbg \
                -a /scratch/mtg_\\\${id}/annotation.row_diff_coord.annodbg \
                ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                2>&1; \
        rm -r /scratch/mtg_\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/megablast;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "query_[1-1]" \
     -o ${DIR}/logs/delta_variants_query_1.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/\\\${id}; \
        cp -H ${DIR}/\\\${id:0:6}/\\\${id}/* /scratch/\\\${id}/; \
        /usr/bin/time -v blastn \
                -query ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -db /scratch/\\\${id}/\\\${id} \
                -use_index true \
                -max_hsps 1 -outfmt '6 qseqid qseq sseq bitscore nident btop' \
                2>&1; \
        rm -r /scratch/\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/blast;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "query_[1-1]" \
     -o ${DIR}/logs/delta_variants_query_1.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=19000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/blast_\\\${id}; \
        cp -H ${DIR}/\\\${id:0:6}/\\\${id}/* /scratch/blast_\\\${id}/; \
        /usr/bin/time -v blastn \
                -query ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -db /scratch/blast_\\\${id}/\\\${id} \
                -max_hsps 1 -outfmt '6 qseqid qseq sseq bitscore nident btop' \
                2>&1; \
        rm -r /scratch/blast_\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/pufferfish_sparse;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "query_[1-1]" \
     -o ${DIR}/logs/delta_variants_query_1.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=200000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/pf_\\\${id}; \
        cp -r ${DIR}/\\\${id:0:6}/\\\${id}/\\\${id} /scratch/pf_\\\${id}; \
        /usr/bin/time -v ~/pufferfish lookup \
                -r ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -i /scratch/pf_\\\${id}/\\\${id}/ \
                2>&1; \
        rm -r /scratch/pf_\\\${id}; \
    done"


DIR=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/pufferfish_sparse;
list=~/metagenome/finished_projects/counting_dbg/hifi_sra/viruses_hifi_data/list.txt;

bsub -J "align_[1-1]" \
     -o ${DIR}/logs/delta_variants_align_1.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=80000] span[hosts=1] select[model==XeonGold_6140]" \
    "for i in \\\$(seq \$((\${LSB_JOBINDEX} * 37 - 37 + 1)) \$((\${LSB_JOBINDEX} * 37))); do \
        file=\\\$(sed -n \\\${i}p $list); \
        id=\\\$(basename \\\${file%.fastq.gz}); \
        mkdir /scratch/pf_align_\\\${id}; \
        cp -r ${DIR}/\\\${id:0:6}/\\\${id}/\\\${id} /scratch/pf_align_\\\${id}; \
        /usr/bin/time -v ~/pufferfish align \
                --read ${DIR}/../mtg/joint/query/redo/delta_variants.fa \
                -i /scratch/pf_align_\\\${id}/\\\${id}/ \
                --genomicReads --primaryAlignment -o /dev/null \
                2>&1; \
        rm -r /scratch/pf_align_\\\${id}; \
    done"

#### Joint index

````bash
K=31;
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data/mtg/joint;
mkdir $DIR;
mkdir $DIR/logs;
list=~/metagenome/data/hifi_sra/viruses_hifi_data/list.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
NUM_THREADS=72;

bsub -J "build_${K}" \
     -o ${DIR}/logs/construct_joint_${K}.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "cat $list | /usr/bin/time -v $METAGRAPH build -v -p $NUM_THREADS \
                -k ${K} \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -o ${DIR}/graph; \

        /usr/bin/time -v $METAGRAPH transform -v -p $NUM_THREADS \
                --state small \
                -o ${DIR}/graph_small \
                ${DIR}/graph.dbg; \
    "


K=31
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data/mtg/joint;
mkdir $DIR/columns;
list=~/metagenome/data/hifi_sra/viruses_hifi_data/list.txt;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "annotate_${K}_[1-1033]" \
     -o ${DIR}/logs/annotate_joint_${K}.lsf \
     -W 24:00 \
     -n 9 -R "rusage[mem=19000] span[hosts=1]" \
    "sed -n \$((\${LSB_JOBINDEX} * 148 - 148 + 1)),\$((\${LSB_JOBINDEX} * 148))p $list \
        | /usr/bin/time -v $METAGRAPH annotate -v \
                --separately -p 3 --threads-each 6 \
                --coordinates \
                --anno-filename \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ${DIR}/columns; \
        "

K=31
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data/mtg/joint;
mkdir $DIR/columns;
list=~/metagenome/data/hifi_sra/viruses_hifi_data/list2.txt;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "annotate_${K}" \
     -o ${DIR}/logs/annotate_joint_${K}_rest.lsf \
     -W 24:00 \
     -n 48 -R "rusage[mem=40000] span[hosts=1]" \
    "cat $list \
        | /usr/bin/time -v $METAGRAPH annotate -v \
                --separately -p 4 --threads-each 24 \
                --coordinates \
                --anno-filename \
                --mem-cap-gb 40 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ${DIR}/columns; \
        "


K=31
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data/mtg/joint;
mkdir ${DIR}/rd_columns;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
NUM_THREADS=36;
bsub -J "rd_coord_${K}_1" \
     -o ${DIR}/logs/rd_joint_${K}.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 0 \
                --mem-cap-gb 400 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ${DIR}/rd_columns/out"

bsub -J "rd_coord_${K}_2" \
     -w "rd_coord_${K}_1" \
     -o ${DIR}/logs/rd_joint_${K}.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 1 \
                --mem-cap-gb 400 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ${DIR}/rd_columns/out"

bsub -J "rd_coord_${K}_3" \
     -w "rd_coord_${K}_2" \
     -o ${DIR}/logs/rd_joint_${K}.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                --anno-type row_diff --coordinates \
                --max-path-length 200 \
                --row-diff-stage 2 \
                --mem-cap-gb 400 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ${DIR}/rd_columns/out"

K=31
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data/mtg/joint;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "rd_brwt_coord_${K}_10M" \
     -o ${DIR}/logs/rd_brwt_coord_joint_${K}.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v -p 72 \
                --anno-type row_diff_brwt_coord \
                --greedy --fast --subsample 10000000 \
                --parallel-nodes 10 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ${DIR}/annotation";

bsub -J "relax_rd_brwt_coord_${K}_10M" \
     -w "rd_brwt_coord_${K}_10M" \
     -o ${DIR}/logs/relax_rd_brwt_coord_joint_${K}.lsf \
     -W 72:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v -p 36 \
                --relax-arity 32 \
                ${DIR}/annotation.row_diff_brwt_coord.annodbg \
                -o ${DIR}/annotation.relaxed";

K=31
DIR=~/metagenome/data/hifi_sra/viruses_hifi_data/mtg/joint;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5/metagraph;
bsub -J "rd_brwt_coord_${K}" \
     -o ${DIR}/logs/rd_brwt_coord_joint_${K}_0.lsf \
     -W 120:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v -p 72 \
                --anno-type row_diff_brwt_coord \
                --parallel-nodes 10 \
                --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                -i ${DIR}/graph.dbg \
                -o ${DIR}/annotation_0";

bsub -J "relax_rd_brwt_coord_${K}" \
     -w "rd_brwt_coord_${K}" \
     -o ${DIR}/logs/relax_rd_brwt_coord_joint_${K}_0.lsf \
     -W 72:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v -p 72 \
                --relax-arity 32 \
                ${DIR}/annotation_0.row_diff_brwt_coord.annodbg \
                -o ${DIR}/annotation_0.relaxed";
```


## Index with k-mer coordinates (lossless read compression)

### Compress with tools for seq compression
```bash
DIR=~/metagenome/data/kingsford;
# mkdir $DIR/compressed;
# mkdir $DIR/compressed/logs;
# find ~/metagenome/raw_data/kingsford/data_fasta -name "*.fasta.gz" > $DIR/compressed/list.txt;
list=$DIR/compressed/list.txt;
bsub -J "gzip[1-$(cat $list | wc -l)]" \
     -o ~/metagenome/data/kingsford/compressed/logs/gzip.lsf \
     -W 24:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\${file%.fasta.gz}); \
          /usr/bin/time -v zcat \\\$file | sed '/^>/d' | tr -d '\n' | gzip -9 > $DIR/compressed/\\\${id}.fasta.gz; \
          /usr/bin/time -v zcat \\\$file | sed '/^>/!d' | tr -d '\n' | gzip -9 > $DIR/compressed/\\\${id}.fasta_headers.gz; \
          /usr/bin/time -v zcat \\\$file | sed 's/^>.*/>/' | awk '/^>/ { if(NR>1) print \\\"\\\"; printf(\\\"%s\n\\\",\\\$0); next; } { printf(\\\"%s\\\",\\\$0);}  END {printf(\\\"\n\\\");}' | gzip -9 > $DIR/compressed/\\\${id}_no_header.fasta.gz"

bsub -J "count_bp[1-$(cat $list | wc -l)]" \
     -w "done(gzip[*])" \
     -o ~/metagenome/data/kingsford/compressed/logs/count_bp.lsf \
     -W 4:00 \
     -n 1 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\${file%.fasta.gz}); \
          /usr/bin/time -v zcat \\\$file | sed '/^>/d' | tr -d '\n' | wc -c > $DIR/compressed/\\\${id}.num_bp;"

bsub -J "spring[1-$(cat ${list} | wc -l)]" \
     -o ~/metagenome/data/kingsford/compressed/logs/noheader_spring.lsf \
     -W 120:00 \
     -n 18 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\${file%.fasta.gz}); \
          rm -rf $DIR/compressed/\\\${id}_no_header_temp; \
          mkdir $DIR/compressed/\\\${id}_no_header_temp; \
          /usr/bin/time -v spring -c -g --fasta-input -t 36 \
                                  -i $DIR/compressed/\\\${id}_no_header.fasta.gz \
                                  -w $DIR/compressed/\\\${id}_no_header_temp \
                                  -o $DIR/compressed/\\\${id}_no_header.spring; \
          rm -r $DIR/compressed/\\\${id}_no_header_temp"

bsub -J "check_spring[1-$(cat ${list} | wc -l)]" \
     -o ~/metagenome/data/kingsford/compressed/logs/check_spring.lsf \
     -W 24:00 \
     -n 4 -R "rusage[mem=10000] span[hosts=1]" \
          "file=\\\$(sed -n \${LSB_JOBINDEX}p ${list}); \
          id=\\\$(basename \\\${file%.fasta.gz}); \
          rm -rf $DIR/compressed/\\\${id}_no_header_temp; \
          mkdir $DIR/compressed/\\\${id}_no_header_temp; \
          rm -rf $DIR/compressed/\\\${id}_no_header.spring.decompressed; \
          /usr/bin/time -v spring -d -t 8 \
                                  -i $DIR/compressed/\\\${id}_no_header.spring \
                                  -w $DIR/compressed/\\\${id}_no_header_temp \
                                  -o $DIR/compressed/\\\${id}_no_header.spring.decompressed; \
          D=\\\$(diff <(zcat $DIR/compressed/\\\${id}_no_header.fasta.gz) $DIR/compressed/\\\${id}_no_header.spring.decompressed | wc -l); \
          echo \\\$D \\\${id} >> ~/metagenome/data/kingsford/compressed/logs/check_spring.out; \
          rm -r $DIR/compressed/\\\${id}_no_header_temp"
```

### with Metagraph
```bash
K=31
DIR=~/metagenome/finished_projects/counting_dbg/kingsford_${K}_coordinates_fork_opt_new;
# git checkout 7a9027fa8c6c29742c7885f77f90d414c39c5b53
rm -rf $DIR;
mkdir $DIR;
mkdir $DIR/logs;

list=~/metagenome/finished_projects/counting_dbg/kingsford_compressed/list.txt;

# rm ~/metagenome/data/kingsford/compressed/list2.txt;
# for file in $(cat ${list}); do
#     if [[ ! -f ${DIR}/$(basename $file)/rd_columns/annotation.column.annodbg.coords ]];
#         then echo $file >> ~/metagenome/data/kingsford/compressed/list2.txt;
#     fi;
# done
# list=~/metagenome/data/kingsford/compressed/list2.txt;

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5_test/metagraph;
NUM_THREADS=8;
bsub -J "build_graph_${K}_[1-$(cat $list | wc -l)]" \
     -o ${DIR}/logs/construct_${K}.lsf \
     -W 24:00 \
     -n 4 -R "rusage[mem=19000] span[hosts=1]" \
NUM_THREADS=36;
bsub -J "build_graph_${K}_rerun_[1-$(cat $list | wc -l)]" \
     -o ${DIR}/logs/construct_${K}_rerun.lsf \
     -w "exit(build_graph_${K}_[*])" \
     -W 24:00 \
     -n 23 -R "rusage[mem=19000] span[hosts=1]" \
NUM_THREADS=36;
bsub -J "build_graph_${K}_rerun2_[1-$(cat $list | wc -l)]" \
     -o ${DIR}/logs/construct_${K}_rerun2.lsf \
     -w "exit(build_graph_${K}_[*]) && exit(build_graph_${K}_rerun_[*])" \
     -W 48:00 \
     -n 36 -R "rusage[mem=25000] span[hosts=1]" \
NUM_THREADS=36;
bsub -J "build_graph_${K}_rerun3_[1-$(cat $list | wc -l)]" \
     -o ${DIR}/logs/construct_${K}_rerun3.lsf \
     -w "exit(build_graph_${K}_[*]) && exit(build_graph_${K}_rerun_[*]) && exit(build_graph_${K}_rerun2_[*])" \
     -W 48:00 \
     -n 36 -R "rusage[mem=40000] span[hosts=1]" \
    "
    file=\\\$(sed -n \${LSB_JOBINDEX}p $list); \
    L=\\\$(zless \\\$file | head -n 1 | grep -Eo '[0-9]+$'); \
    if [[ \\\$L -le 30 ]]; then exit 0; fi; \
    id=\\\$(basename \\\$file); \
    rm -r ${DIR}/\\\$id; \
    mkdir ${DIR}/\\\$id; \
    mkdir ${DIR}/\\\$id/logs; \

    /usr/bin/time -v $METAGRAPH build -v -p $NUM_THREADS \
            -k ${K} \
            --mem-cap-gb 40 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -o ${DIR}/\\\${id}/graph \
            \\\$file \
            2> ${DIR}/\\\${id}/logs/build.log; \

    /usr/bin/time -v $METAGRAPH transform -v -p $NUM_THREADS \
            --state small \
            --index-ranges 1 \
            -o ${DIR}/\\\${id}/graph_small \
            ${DIR}/\\\${id}/graph.dbg \
            2> ${DIR}/\\\${id}/logs/transform.log; \

    /usr/bin/time -v $METAGRAPH annotate -v -p $NUM_THREADS \
            --coordinates \
            --anno-filename \
            --mem-cap-gb 40 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/\\\${id}/graph.dbg \
            -o ${DIR}/\\\${id}/annotation \
            \\\$file \
            2> ${DIR}/\\\${id}/logs/annotate.log; \

    mkdir ${DIR}/\\\${id}/rd_columns; \

    /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
            --anno-type row_diff --coordinates \
            --max-path-length 200 \
            --row-diff-stage 0 \
            --mem-cap-gb 200 \
            --num-kmers-in-seq \\\$((\\\${L}-${K}+1)) \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/\\\${id}/graph.dbg \
            -o ${DIR}/\\\${id}/rd_columns/out \
            ${DIR}/\\\${id}/annotation.column.annodbg \
            2> ${DIR}/\\\${id}/logs/coord_rd_1.log; \

    /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
            --anno-type row_diff --coordinates \
            --max-path-length 200 \
            --row-diff-stage 1 \
            --mem-cap-gb 200 \
            --num-kmers-in-seq \\\$((\\\${L}-${K}+1)) \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/\\\${id}/graph.dbg \
            -o ${DIR}/\\\${id}/rd_columns/out \
            ${DIR}/\\\${id}/annotation.column.annodbg \
            2> ${DIR}/\\\${id}/logs/coord_rd_2.log; \

    /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
            --anno-type row_diff --coordinates \
            --max-path-length 200 \
            --row-diff-stage 2 \
            --mem-cap-gb 200 \
            --num-kmers-in-seq \\\$((\\\${L}-${K}+1)) \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/\\\${id}/graph.dbg \
            -o ${DIR}/\\\${id}/rd_columns/out \
            ${DIR}/\\\${id}/annotation.column.annodbg \
            2> ${DIR}/\\\${id}/logs/coord_rd_3.log; \

    rm ${DIR}/\\\${id}/annotation.column.annodbg; \
    rm ${DIR}/\\\${id}/annotation.column.annodbg.coords; \
    rm ${DIR}/\\\${id}/rd_columns/annotation.row_count; \
    rm ${DIR}/\\\${id}/rd_columns/annotation.row_reduction; \
    rm ${DIR}/\\\${id}/graph.dbg.succ; \
    rm ${DIR}/\\\${id}/graph.dbg.succ_boundary; \
    rm ${DIR}/\\\${id}/graph.dbg.pred; \
    rm ${DIR}/\\\${id}/graph.dbg.pred_boundary"; \


for file in '~/metagenome/data/coordinates_K/data/SRR13577847_subreads.fastq.gz' ; do
    bsub -J "gzip" \
         -o /dev/null \
         -W 24:00 \
         -n 1 -R "rusage[mem=10000] span[hosts=1]" \
              "file=${file}; \
              id=\\\$(basename \\\${file%.gz}); \
              /usr/bin/time -v zcat ${file} | sed -n '1~4s/^@/>/p;2~4p' | sed 's/^>.*/>/' | gzip -9 > ~/metagenome/data/coordinates_K/data/\\\${id}_no_header.fasta.gz"
done

for file in '~/metagenome/finished_projects/counting_dbg/coordinates_K/data/SRR11304401_subreads.fastq.gz' \
            '~/metagenome/finished_projects/counting_dbg/coordinates_K/data/SRR13684276.fastq.gz' \
            '~/metagenome/finished_projects/counting_dbg/coordinates_K/data/SRR13577847_subreads.fastq.gz' \
            '/cluster/work/grlab/projects/metagenome/data/alignment/completeness/wgs_samples/fq/fq2/SRR4063132/SRR4063132_subreads.fastq.gz' \
            '/cluster/work/grlab/projects/metagenome/data/alignment/completeness/wgs_samples/fq/fq2/SRR386922/SRR386922_subreads.fastq.gz' \
            '/cluster/work/grlab/projects/metagenome/data/alignment/completeness/wgs_samples/fq/fq2/SRR3747284/SRR3747284_subreads.fastq.gz' \
            '/cluster/work/grlab/projects/metagenome/data/alignment/completeness/wgs_samples/fq/fq2/SRR4235456/SRR4235456_subreads.fastq.gz' \
            '/cluster/work/grlab/projects/metagenome/data/alignment/completeness/wgs_samples/fq/fq2/SRR3747411/SRR3747411_subreads.fastq.gz' \
            '~/metagenome/finished_projects/counting_dbg/kingsford_compressed/SRR805801_no_header.fasta.gz' \
            '/cluster/work/grlab/projects/metagenome/raw_data/human/HG002/PacBio_SequelII_CCS_11kb/m64011_181218_235052.fastq.gz' ; do
    for K in {31,}; do
        DIR=~/metagenome/data/coordinates_K/$K;
        # rm -rf $DIR;
        mkdir -p $DIR;
        mkdir -p $DIR/logs;

        METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_dna5_test/metagraph;
        NUM_THREADS=16;
        # bsub -J "construct_${K}_[1-10]" \
        bsub -J "build_graph_${K}" \
             -o ${DIR}/logs/construct_${K}.lsf \
             -W 24:00 \
             -n 8 -R "rusage[mem=40000] span[hosts=1]" \
            "
            file=$file; \
            id=\\\$(basename \\\$file); \
            mkdir ${DIR}/\\\$id; \
            mkdir ${DIR}/\\\$id/logs; \

            /usr/bin/time -v $METAGRAPH build -v -p $NUM_THREADS \
                    -k ${K} \
                    --mem-cap-gb 40 \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -o ${DIR}/\\\${id}/graph \
                    \\\$file; \

            /usr/bin/time -v $METAGRAPH transform -v -p $NUM_THREADS \
                    --state small \
                    --index-ranges 3 \
                    -o ${DIR}/\\\${id}/graph_small \
                    ${DIR}/\\\${id}/graph.dbg; \

            /usr/bin/time -v $METAGRAPH annotate -v -p $NUM_THREADS \
                    --coordinates \
                    --anno-filename \
                    -i ${DIR}/\\\${id}/graph.dbg \
                    -o ${DIR}/\\\${id}/annotation \
                    \\\$file; \

            mkdir ${DIR}/\\\${id}/rd_columns; \

            /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                    --anno-type row_diff --coordinates \
                    --max-path-length 200 \
                    --row-diff-stage 0 \
                    --mem-cap-gb 200 \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -i ${DIR}/\\\${id}/graph.dbg \
                    -o ${DIR}/\\\${id}/rd_columns/out \
                    ${DIR}/\\\${id}/annotation.column.annodbg; \

            /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                    --anno-type row_diff --coordinates \
                    --max-path-length 200 \
                    --row-diff-stage 1 \
                    --mem-cap-gb 200 \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -i ${DIR}/\\\${id}/graph.dbg \
                    -o ${DIR}/\\\${id}/rd_columns/out \
                    ${DIR}/\\\${id}/annotation.column.annodbg; \

            /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
                    --anno-type row_diff --coordinates \
                    --max-path-length 200 \
                    --row-diff-stage 2 \
                    --mem-cap-gb 200 \
                    --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
                    -i ${DIR}/\\\${id}/graph.dbg \
                    -o ${DIR}/\\\${id}/rd_columns/out \
                    ${DIR}/\\\${id}/annotation.column.annodbg; \

            rm ${DIR}/\\\${id}/annotation.column.annodbg; \
            rm ${DIR}/\\\${id}/annotation.column.annodbg.coords; \
            rm ${DIR}/\\\${id}/rd_columns/annotation.row_count; \
            rm ${DIR}/\\\${id}/rd_columns/annotation.row_reduction; \
            rm ${DIR}/\\\${id}/graph.dbg.succ; \
            rm ${DIR}/\\\${id}/graph.dbg.succ_boundary; \
            rm ${DIR}/\\\${id}/graph.dbg.pred; \
            rm ${DIR}/\\\${id}/graph.dbg.pred_boundary"; \
    done
done
```



```bash
list=~/metagenome/data/row_diff/subsets/refseq/8.txt
DIR=~/metagenome/finished_projects/counting_dbg/refseq_fungi_coord;
mkdir $DIR;
mkdir $DIR/logs;

cat ~/metagenome/data/row_diff/subsets/refseq/8.txt | xargs -P 100 -I {} sh -c "zcat {} | gzip -9 > {}.9"
cat ~/metagenome/data/row_diff/subsets/refseq/8.txt | xargs -P 200 -I {} sh -c "zcat {} | grep '>' >> ~/metagenome/data/refseq_fungi_coord/headers.txt"
cat ~/metagenome/data/row_diff/subsets/refseq/8.txt | xargs -P 200 -I {} sh -c "zcat {} | grep -v '>' | tr -d '\n' | wc -c" | awk "{sum+=\$1}END{print sum}"
cat ~/metagenome/data/row_diff/subsets/refseq/8.txt | xargs -I {} sh -c "ls -l {}.9" | sizeb
cat ~/metagenome/data/refseq_fungi_coord/logs/annotate_graph.lsf | grep "Number of coordinates" | cut -d' ' -f7 | awk "{sum+=\$1}END{print sum}"

DIR=~/metagenome/data/refseq_fungi_coord;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
bsub -J "fungi_graph" \
     -oo $DIR/logs/build_graph.lsf \
     -W 12:00 \
     -n 15 -R "rusage[mem=19000] span[hosts=1]" \
    "cat ${list} \
        | /usr/bin/time -v $METAGRAPH build -v \
            -k 31 \
            -p 30 \
            --mem-cap-gb 300 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -o $DIR/graph"; \

mkdir $DIR/columns;
bsub -J "annotate_fungi" \
     -w "fungi_graph" \
     -oo ${DIR}/logs/annotate_graph.lsf \
     -W 12:00 \
     -n 18 -R "rusage[mem=15000] span[hosts=1]" \
    "cat ${list} \
        | /usr/bin/time -v $METAGRAPH annotate \
            -i $DIR/graph.dbg \
            --anno-header \
            --separately \
            --coordinates \
            -o ${DIR}/columns \
            -p 36"; \

DIR=~/metagenome/data/refseq_fungi_coord;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_master/metagraph;

mkdir $DIR/rd;
mkdir $DIR/rd/rd_columns;
ln -s $DIR/graph.dbg ${DIR}/rd/graph.dbg;

DIR=~/metagenome/data/refseq_fungi_coord;
bsub -J "fungi_rd_0" \
     -oo ${DIR}/logs/rd_coord_0.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --coordinates \
            --row-diff-stage 0 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72";

DIR=~/metagenome/data/refseq_fungi_coord;
bsub -J "fungi_rd_1" \
     -w "fungi_rd_0" \
     -oo ${DIR}/logs/rd_coord_1.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --coordinates \
            --row-diff-stage 1 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72";

DIR=~/metagenome/data/refseq_fungi_coord;
bsub -J "fungi_rd_2" \
     -w "fungi_rd_1" \
     -oo ${DIR}/logs/rd_coord_2.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff \
            --coordinates \
            --row-diff-stage 2 \
            --mem-cap-gb 500 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/rd/rd_columns/out \
            -p 72";

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/refseq_fungi_coord;
bsub -J "fungi_rd_coord" \
     -w "fungi_rd_2" \
     -oo ${DIR}/logs/rd_coord.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=4000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_coord \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72";

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_test/metagraph;
DIR=~/metagenome/data/refseq_fungi_coord;
bsub -J "fungi_rd_brwt_coord" \
     -oo ${DIR}/logs/rd_brwt_coord.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_brwt_coord \
            --greedy --fast --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72";

METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
DIR=~/metagenome/data/refseq_fungi_coord;
bsub -J "fungi_rd_brwt_coord_relax" \
     -oo ${DIR}/logs/rd_brwt_coord_relax.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH relax_brwt -v \
            --relax-arity 8096 \
            -o ${DIR}/annotation.relaxed \
            ${DIR}/annotation.row_diff_brwt_coord.annodbg \
            -p 72";





K=21
DIR=~/metagenome/data/coordinates_K/data/random_counts;
#rm -rf $DIR;
mkdir $DIR;
mkdir $DIR/logs;
METAGRAPH=~/projects/projects2014-metagenome/metagraph/build_release/metagraph;
NUM_THREADS=36;
bsub -J "build_graph" \
     -oo ${DIR}/logs/construct.lsf \
     -W 24:00 \
     -n 36 -R "rusage[mem=19000] span[hosts=1]" \
    "/usr/bin/time -v $METAGRAPH build -v -p $NUM_THREADS \
            -k ${K} \
            --mem-cap-gb 40 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -o ${DIR}/graph \
            ${DIR}/random_1M.fasta.gz \
            2> ${DIR}/logs/build.log; \

    /usr/bin/time -v $METAGRAPH transform -v -p $NUM_THREADS \
            --state small \
            --index-ranges 1 \
            -o ${DIR}/graph_small \
            ${DIR}/graph.dbg \
            2> ${DIR}/logs/transform.log; \

    /usr/bin/time -v $METAGRAPH annotate -v -p $NUM_THREADS \
            --count-kmers \
            --count-width 32 \
            --anno-filename \
            --mem-cap-gb 40 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/graph.dbg \
            -o ${DIR}/annotation \
            ${DIR}/random_1M.fasta.gz \
            2> ${DIR}/logs/annotate.log; \

    /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type int_brwt \
            --greedy --fast --subsample 1000000 \
            -o ${DIR}/annotation \
            ${DIR}/annotation.column.annodbg \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/count_brwt.log; \

    mkdir ${DIR}/rd_columns; \

    /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
            --anno-type row_diff --count-kmers \
            --max-path-length 200 \
            --row-diff-stage 0 \
            --mem-cap-gb 200 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/graph.dbg \
            -o ${DIR}/rd_columns/out \
            ${DIR}/annotation.column.annodbg \
            2> ${DIR}/logs/coord_rd_1.log; \

    /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
            --anno-type row_diff --count-kmers \
            --max-path-length 200 \
            --row-diff-stage 1 \
            --mem-cap-gb 200 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/graph.dbg \
            -o ${DIR}/rd_columns/out \
            ${DIR}/annotation.column.annodbg \
            2> ${DIR}/logs/coord_rd_2.log; \

    /usr/bin/time -v $METAGRAPH transform_anno -v -p $NUM_THREADS \
            --anno-type row_diff --count-kmers \
            --max-path-length 200 \
            --row-diff-stage 2 \
            --mem-cap-gb 200 \
            --disk-swap ~/metagenome/scratch/nobackup/stripe_1 \
            -i ${DIR}/graph.dbg \
            -o ${DIR}/rd_columns/out \
            ${DIR}/annotation.column.annodbg \
            2> ${DIR}/logs/coord_rd_3.log; \

    rm ${DIR}/rd_columns/annotation.row_count; \
    rm ${DIR}/rd_columns/annotation.row_reduction; \
    rm ${DIR}/graph.dbg.succ; \
    rm ${DIR}/graph.dbg.succ_boundary; \
    rm ${DIR}/graph.dbg.pred; \
    rm ${DIR}/graph.dbg.pred_boundary; \

    find ${DIR}/rd/rd_columns -name \"*.column.annodbg\" \
        | /usr/bin/time -v $METAGRAPH transform_anno -v \
            --anno-type row_diff_int_brwt \
            --greedy --fast --subsample 1000000 \
            -i ${DIR}/rd/graph.dbg \
            -o ${DIR}/annotation \
            -p 72 --parallel-nodes 10 \
            2>&1 | tee ${DIR}/logs/count_rd_brwt.log"; \
```
