# Data analysis

The following is a description of the data analysis steps performed on the raw
sequencing output. It aims to provide all necessary information to go from the
raw FASTQ files output by the Illumina sequencing machines to the processed data
used to build the graphs.

Some tasks were performed on the NYU HPC platform (running the
 [slurm](https://slurm.schedmd.com/) workload manager).

### Included analyses

* [Downloading raw data](#downloading-raw-data)
* [Read quality trimming](#read-quality-trimming)
* [Preparation of hybrid SK1-S288C yeast genome](README.md#preparation-of-hybrid-s1-s288c-yeast-genome)
* [List of SNPs between SK1 and S288C](README.md#list-of-snps-between-sk1-and-s288c)
* [Read alignment, pileup and peak calling pipeline](README.md#read-alignment-pileup-and-peak-calling-pipeline)
* [Summary of aligned read coverage per genomic position](README.md#summary-of-aligned-read-coverage-per-genomic-position)

## Downloading raw data

The raw sequencing data can be obtained from NCBI's Gene Expression Omnibus through GEO Series accession number [GSE114731](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114731) (private until publication). The following code uses NCBI's Entrez Direct and SRA Toolkit (see [Download Tools](https://www.ncbi.nlm.nih.gov/home/tools/)) to download all SRR files in the project and dump them to FASTQ files (will work after publication of the paper, once the submission is public).

> _Please note that the follownig code will download several relatively large files. Alternatively, you can download selected files by using their individual SRA code. For example, download "SK1 WT / S288c WT Input rep 1" by running:_ `$ fasterq-dump SRX4108751`

```bash
# Use GSE114731's BioProject name in query
esearch -db sra -query PRJNA472390 | \
efetch --format runinfo | cut -d ',' -f 1 | \
grep SRR | xargs fasterq-dump
```

## Read quality trimming

Read quality trimming was performed with
[Trim Galore!](https://github.com/FelixKrueger/TrimGalore) according to the
following example.


```bash
srun -t3:00:00 --mem=4000 --pty /bin/bash
module load cutadapt/intel/1.12
module load fastqc/0.11.5
module load trim_galore/0.4.4

mkdir Trim-galore_output

for FILE in *.fastq.gz
do
    echo
    echo "--------------"
    echo ">>>>> Trimming ${FILE}"
    trim_galore -q 28 --fastqc --output_dir Trim-galore_output ${FILE}
done
```


## Preparation of hybrid SK1-S288C yeast genome

The hybrid genome is prepared by simply concatenating both reference genomes
end-to-end, to generate a new combined reference genome. The following code is
the code used for this task.

```bash
mkdir S288C_SK1_Yue_hybrid_genome && cd S288C_SK1_Yue_hybrid_genome

# Download nuclear genomes
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/SK1.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/SK1.all_feature.gff.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/S288c.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/S288c.all_feature.gff.gz
gunzip *

# SK1 FASTA has 32 lines, S288c has 221656;
# Unlikely to make any difference, but delete excess new lines
# in S288c.genome.fa for consistency
awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n }' S288c.genome.fa \
>> new_S288c.genome.fa
mv new_S288c.genome.fa S288c.genome.fa

# Add strain name to chr names (to make them distinguishable)
sed -i -E 's/(>chr[IVX]+)/\1_S288C/' S288c.genome.fa
sed -i -E 's/(>chr[IVX]+)/\1_SK1/' SK1.genome.fa

# Concatenate genomes (cat alone would not introduce new line between files)
cat S288c.genome.fa <(echo) SK1.genome.fa > S288c_SK1_Yue.fa
```


## List of SNPs between SK1 and S288C

Align SK1 and S288C genome assemblies
([Yue _et al._, Nat Genet 2017](https://www.ncbi.nlm.nih.gov/pubmed/28416820))
and get SNPs using [Mauve](http://darlinglab.org/mauve/mauve.html)
([Darling _et al._, Genome Res 2014](https://www.ncbi.nlm.nih.gov/pubmed/15231754))
as described below.

#### Download genomes

```bash
mkdir Yue_SK1_v_S288C_SNPs
cd Yue_SK1_v_S288C_SNPs

wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/SK1.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/SK1.all_feature.gff.gz

wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_Genome/S288c.genome.fa.gz
wget http://yjx1217.github.io/Yeast_PacBio_2016/data/Nuclear_GFF/S288c.all_feature.gff.gz

gunzip *
```

#### Generate list of SNPs

* Select **`Align with progressiveMauve...`** in the Mauve GUI
    * Select genome files and define output file name:
        - `S288c.genome.fa`
        - `SK1.genome.fa`
        - Output file name: `S288c_v_SK1.aln`
* Select **`Tools > Export > Export SNPs`** in the Mauve GUI
    * Define output file name:
        - Output file names: `S288c_v_SK1.snp`

#### Format SNP table for use

Run the following R code.

```r
# Load SNPs
S288cvSK1_snp <- readr::read_tsv(here('data/S288c_v_SK1.snp'))
head(S288cvSK1_snp, 10)

message('Number of identified SNPs: ', nrow(S288cvSK1_snp))

# SNPs on different contigs:
S288cvSK1_snp_diff_contig <- S288cvSK1_snp[S288cvSK1_snp$sequence_1_Contig != S288cvSK1_snp$sequence_2_Contig, ]

head(S288cvSK1_snp_diff_contig)

message('Number of SNPs on different chromosomes: ', nrow(S288cvSK1_snp_diff_contig))
message('(Remaining SNPs on matching chromosomes: ',
        nrow(S288cvSK1_snp[S288cvSK1_snp$sequence_1_Contig == S288cvSK1_snp$sequence_2_Contig, ]))

# Drop the 236 non-matched SNPs
S288cvSK1_snp <- S288cvSK1_snp[S288cvSK1_snp$sequence_1_Contig == S288cvSK1_snp$sequence_2_Contig, ]

# Make tidy table
SNPs <- S288cvSK1_snp[, c('sequence_1_Contig', 'sequence_1_PosInContg',
                          'sequence_2_Contig', 'sequence_2_PosInContg')]

SNPs$sequence_1_Contig <- paste0(SNPs$sequence_1_Contig, '_S288C')
SNPs$sequence_2_Contig <- paste0(SNPs$sequence_2_Contig, '_SK1')

colnames(SNPs) <- c('chr_S288C', 'position_S288C', 'chr_SK1', 'position_SK1')

S288C <- subset(SNPs, select=c(chr_S288C, position_S288C))
SK1 <- subset(SNPs, select=c(chr_SK1, position_SK1))
colnames(S288C) <- c('chr', 'position')
colnames(SK1) <- c('chr', 'position')
SNPs <- rbind(S288C, SK1)

head(SNPs)

# Save new table to file
readr::write_tsv(SNPs, here('data/S288c_v_SK1.snp'))
```


## Read alignment, pileup and peak calling pipeline

Ran pipeline using the `slurm` file `ChIPseq_Pipeline_hybrid_genome.sbatch` (see
[ChIP-seq Pipeline](ChIPseq_Pipeline_hybrid_genome/README.md)
section). The
following is the submission code template.

```bash
sbatch --export EXPID="SampleID_SK1_S288C_Yue_PM_SPMR",\
RUNDIR="Run-dir/",\
CHIP="Run-dir/chip_trimmed.fastq",\
INPUT="Run-dir/input_trimmed.fastq",\
GENNAME="Run-dir/S288C_SK1_Yue_hybrid_genome/S288c_SK1_Yue" \
ChIPseq_Pipeline_hybrid_genome.sbatch
```


## Summary of aligned read coverage per genomic position

Read coverage at each genomic position (or read pileup) was obtained using
`bedtools genomecov`. It was used, for example, for an alternative way of
computing spike-in normalization factors, to compare with the standard way
(using total number of mapped reads).

The following is example code that can be run in a directory containing aligned
read maps in BAM format (sorted by position) to generate coverage files in
BedGraph format.

```bash
for FILE in *_sorted.bam
do
    bedtools genomecov -ibam ${FILE} -bg > Pileup_${FILE%_sorted.bam}.bdg
done
```

The generated bedtools pileups were used to calculate spike-in normalization
factors using different input data. The R code can be found in file
`helper_spikein_normalization_factor_alternative_calculations.R`. The final
results table was saved to file
 `data/spikein_normalization_factors_using-different_methods.csv`.
