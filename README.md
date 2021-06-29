# genomicSeq_callSNPs

simple folder to teach Courtney and Risa how to call SNPs in genomic sequence (Drosophila)

## github location (notes for Janet)

public:
https://github.com/jayoung/genomicSeq_callSNPs

local:
/fh/fast/malik_h/user/jayoung/forOtherPeople/forCourtney/genomicSeq_callSNPs

(a reminder of how to sync to github: `git add --all .` then `git commit` then `git push`)

## download reads to Hutch server

Done.

For each DNA sample, Novogene gave us two files, containing the forward and reverse reads (sample1_1.fastq.gz and sample1_2.fastq.gz). 

Some other sequencing facilities (like the Hutch facility) will give you multiple files for each sample+direction combination, so you have to decide when to combine the files that represent the same sample+direction. Many programs can take multiple input files and consider them as part of the same dataset.


## run fastqc (quality control, get some basic stats on the data)

For each file of reads (e.g. sample1_1.fastq.gz), we want to look at quality control metrics.

On the command-line, run a command that looks like this - you will change the last two parts for each sample+direction read dataset. You'll want to specify a uniquely named output directory (--outdir) as well as the input fastq.gz file.

The general format of the fastqc command is `fastqc [options] input file(s)` but I give you more specifics below.

```
module load FastQC/0.11.9-Java-11 

fastqc --format fastq --threads 4 --contaminants /fh/fast/malik_h/user/jayoung/general_notes/NextGenSeqStuff/adapterSequences/variousAdaptersBothStrands.fa.txt --outdir fastqc_output/Sample_1_Control_DNA/R2 ../FASTQ_files/Sample_1_Control_DNA/mySample1_1.fq.gz

module purge
```

it is probably self-evident what each of the options mean, but if you want to learn more, you can (after doing the module load), do `fastqc --help`  (or google, of course)

## maybe make some very small fastq files so we can test whether our code works as expected.

These commands would take the first 400 lines of the forward and reverse files and save as new files. Each sequence read takes 4 lines, so 400 lines = 100 sequences.

```
zcat mySample1_1.fq.gz | head -400 > mySample1_1.first100.fq
gzip mySample1_1.first100.fq

zcat mySample1_2.fq.gz | head -400 > mySample1_2.first100.fq
gzip mySample1_2.first100.fq
```

## trim reads, maybe (depends on what fastqc output looks like)

Remove low quality regions, adapter sequences

Input and output files are both in fastq.gz format

The cutadapt program can probably do what we want it to.  We'll run it on BOTH forward and reverse read files at the same time.  We want the two files to stay synced up, i.e. to always contain paired reads in the same order, so if cutadapt decides to reject the forward read, it'll be able to reject the reverse read as well.

The command will look roughly like this, but we'll need to modify it depending on how fastqc output looks. We would replace ADAPT1 and ADAPT2 with the sequences of adapters we want to trim.

```
module load cutadapt/2.9-foss-2019b-Python-3.7.4 

# to check out all the options:
cutadapt --help

# our actual command will look something like this
cutadapt [any_more_options] --quality-cutoff 30 --cores 4 -a ADAPT1 -A ADAPT2 -o mySample1_1.trimmed.fq.gz -p mySample1_2.trimmed.fq.gz mySample1_1.fq.gz mySample1_2.fq.gz > mySample1_2.fq.cutadapt.log.txt 

module purge
```

## map to reference genome

input:   fastq.gz files and reference genome assembly
output:  bam files

For the reference genome, I suggest we use the dm6 version of the Drosophila melanogaster genome assembly.  

I also suggest we use a version of dm6 where I have removed the 'chrUn' and 'random' sequences (chinks of sequence that were not assigned to any chromosome).  Some of those are alternative alleles (I think?) of other regions that are in the 'regular' chromosome sequences, and having >1 version of a chromosome sequence present will confuse bwa, as for reads that map to multiple locations, it chooses one location at random.

To map reads we'll use the BWA program. See http://bio-bwa.sourceforge.net. There are three algorithms in the BWA package - we'll use BWA-MEM, which is the best one for reads >70bp

I already formatted that genome assembly for BWA mapping (using the `bwa index` algorithm): `/fh/fast/malik_h/grp/public_databases/UCSC/fly_Aug2014_dm6/dm6_withoutChrUnRandom/dm6_withoutChrUnRandom.fa_bwaFormat/dm6_withoutChrUnRandom.fa`

The general form of a bwa mem command is `bwa mem [options] refGenome.fa read1.fq read2.fq > alignedReads.sam`

Option -t specifies number of threads (CPUs)

```
module load BWA/0.7.17-GCC-10.2.0

bwa mem -t 4 /fh/fast/malik_h/grp/public_databases/UCSC/fly_Aug2014_dm6/dm6_withoutChrUnRandom/dm6_withoutChrUnRandom.fa_bwaFormat/dm6_withoutChrUnRandom.fa sample1_1.trimmed.fastq.gz sample1_2.trimmed.fastq.gz > sample1.dm6.sam

module purge
```

BWA creates a 'sam' format file that is NOT sorted by genomic position. You can look at this using a command like 'more' or 'less' or 'cat'.   Most downstream programs want a sorted file, and prefer 'bam' format (similar to sam but compressed, so files are smaller), so we will convert to bam and sort, using samtools (see http://www.htslib.org/doc/samtools.html).  We then create a bam index file that helps other programs use the bam file. We'll remove the intermediate files to save disk space.

The -@ option for samtools specifies number of threads.

We'll also run `samtools flagstat` on the bam file to see what % of reads mapped to the genome

```
module load SAMtools/1.11-GCC-10.2.0

samtools view -@ 4 -Sb sample1.dm6.sam > sample1.dm6.bam
samtools sort -@ 4 -O bam sample1.dm6.bam > sample1.dm6.sorted.bam
samtools index sample1.dm6.sorted.bam
samtools flagstat sample1.dm6.sorted.bam > sample1.dm6.sorted.bam.flagstats

rm sample1.dm6.sam  sample1.dm6.bam

module purge
```

I suggest putting the BWA and samtools commands in a single shell script file, one for each of your two samples. 


## use IGV to browse reads in the region of the genes you knocked out

to get ready for this, download and install IGV on your mac:  

https://software.broadinstitute.org/software/igv/download

(I think the 'IGV MacOS App, Java included' version)

We'll start up IGV, load the reference genome and the bam files of mapped reads.

We can talk through how to look at it together, but you might also want to browser the user guide: https://software.broadinstitute.org/software/igv/userguide


### Loading the reference genome, and annotations (genes, repeats, etc)

Option 1: (simplest, but I'm not sure how it will behave)
Use menu option 'Genomes-Load genome from server' and find D. melanogaster (dm6).


Option 2:
Load the reference genome by pulling up menu option 'Genomes-Load genome from file', then navigating to a version of the reference genome that I have already indexed this for use in IGV: `/fh/fast/malik_h/grp/public_databases/UCSC/fly_Aug2014_dm6/dm6_withoutChrUnRandom/dm6_withoutChrUnRandom.fa` 

We'll also load up the gene annotations - now we use the 'File-Load from file' option and navigate to this file (again - I've already indexed this for use in IGV):
`/fh/fast/malik_h/grp/public_databases/UCSC/fly_Aug2014_dm6/misc_tracks/dm6_refGene_2018_Jul9.sorted.changeNames_withoutChrUnRandom.bed`

Perhaps also repetitive element annotations - `/fh/fast/malik_h/grp/public_databases/UCSC/fly_Aug2014_dm6/mask_tracks/dm6_rmsk_2018_Jul9_withoutChrUnRandom.bed`


### Loading your mapped reads (and later, your SNP calls)

Use menu option 'File-Load from file' and find the bam files you made using BWA.

Then navigate to your region of interest. Sometimes it's helpful to use the UCSC genome browser to look up where in the dm6 assembly something is (e.g. Arp53D is at chr2R:16,774,246-16,775,676).  You can copy-paste coordinates into the IGV box to the left of the 'Go' button.

We should be able to see the deleted region around your gene of interest, as well as some reads where the mate-pair did not map to the reference genome

We will also be able to load the vcf files you will make later when you call SNPs




## call SNPs and small indels

'SNP' in this context means difference from the reference genome assembly

we'll use the GATK pipeline - see the 'cohort data' section of this page: 
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-


input: bam files

output: vcf file (variant calling format)


### filter SNPs

ignore low confidence SNPs

ignore SNPs where knockout and paired control are identical (there will be MANY!)

perhaps annotate SNPs with information about which gene they are in, what effect they have on the gene (e.g. synonymous/non-synonymous), what other strains they are found in, 


## call larger variants 

(insertions, deletions, copy number variation, translocations, inversions)

Various other algorithms can do this - we'll want to run several, probably.
