# Strand and be Counted
## Aim
Sometimes it can be tricky to know which strand options to use in RNA-seq processing tools, even when consulting the documentation. *strand_and_be_counted* takes a pair of FASTQ.gz reads, the reference and GTF annotation file, and reports on which strand and orientation options to use for HISAT2 and htseq-count.

## Requirements
* samtools
* htseq
* hisat2
* pandas

These can be installed with conda or mamba, e.g.

```mamba create -n strand_and_be_counted -c bioconda htseq hisat2 pandas samtools```

## Usage

### Example
```
./strand_and_be_counted.py \
    -r ref.fa \
    -g ref.gtf \
    --r1 sample.R1.fastq.gz \
    --r2 sample.R2.fastq.gz \
    --num_reads 100000 \
    --genome_prop 0.1 \
    --out output_dir
```

### Option Descriptions
```
options:
  -h, --help            show this help message and exit
  -r REF, --ref REF     reference FASTA
  -g GTF, --gtf GTF     reference GTF
  --r1 R1               R1.fastq.gz
  --r2 R2               R2.fastq.gz
  --num_reads NUM_READS
                        Number of reads to check. -1 = all [-1]
  --genome_prop GENOME_PROP
                        Proportion of genome to use [1.0]
  --idattr IDATTR       htseq count gene name attribute [gene_id]
  --out OUT             Output directory [./out]
  --threads THREADS     Threads [64]
```

## Analysis steps
*strand_and_be_counted* carries out the following steps
1. Subsamples the read pairs
2. Subsamples the reference genome and GTF
3. Indexes the genome and align reads with HISAT2 using ```--fr/--rf/--ff``` options
4. Selects the best results set and obtains gene counts (exon) with htseq-count, using ```--strand no/reverse/yes``` options
5. Selects the best results set and reports recommended parameters

## Example Output
A file called **results.txt** is generated after each run, which looks like this:
```
*** HISAT2 results
    total_map_rate  concordant_map_rate  proportion_concordant
fr            7.55                 6.86               0.908609
rf            7.55                 1.42               0.188079
ff            7.54                 0.01               0.001326
*** htseq-count results using --fr bam data
         __no_feature  __ambiguous  __too_low_aQual  __not_aligned  __alignment_not_unique  assigned_reads
yes              7579            0                0              0                     427              14
reverse          5038            1                0              0                     427            2554
no               5025            1                0              0                     427            2567

Recommend using --fr for HISAT2
Recommend using -s reverse for htseq-count
```

