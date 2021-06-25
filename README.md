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

We'll use dm6 version of the reference genome

We'll use the BWA program


input:   fastq.gz files

output:  bam files

then we'll run `samtools flagstat` on the bam file to see what % of reads mapped to the genome


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
