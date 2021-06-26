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


## trim reads, maybe (depends on what fastqc output looks like)

remove low quality regions

remove adapter sequences

input and output files are both in fastq.gz format


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




## use IGV to browse reads in the region of the genes you knocked out

to get ready for this, download and install IGV on your mac:  

https://software.broadinstitute.org/software/igv/download

(I think the 'IGV MacOS App, Java included' version)


We'll start up IGV, load the reference genome and the bam files of mapped reads.

We should be able to see the deleted region in your gene of interest, as well as some reads where the mate-pair did not map to the reference genome

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
