# genomicSeq_callSNPs

simple folder to teach Courtney and Risa how to call SNPs in genomic sequence (Drosophila)

https://github.com/jayoung/genomicSeq_callSNPs


## download reads to Hutch server

done


## run fastqc (quality control, get some basic stats on the data)

On the command-line:

```
module load FastQC/0.11.9-Java-11 
fastqc --outdir fastqc_output/Sample_1_Control_DNA/R2 --format fastq --threads 4 --contaminants /fh/fast/malik_h/user/jayoung/general_notes/NextGenSeqStuff/adapterSequences/variousAdaptersBothStrands.fa.txt ../FASTQ_files/Sample_1_Control_DNA/mySample1_1.fq.gz
module purge
```

it is probably self-evident what each of the options mean, but if you want to learn more, you can (after doing the module load), do `fastqc --help`  (or google, of course)


## trim reads, maybe (depends on what fastqc output looks like)

remove low quality regions
remove adapter sequences


## map to reference genome

We'll use dm6 version of the reference genome
We'll use the BWA program

input:   fastq.gz files
output:  bam files

then we'll run `samtools flagstat` on the bam file to see what % of reads mapped to the genome


## use IGV to browse reads in the region of the genes you knocked out


## call SNPs and small indels

GATK

see here, the 'cohort data' section
https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-


## call larger variants 

(insertions, deletions, copy number variation, translocations, inversions)

Various other algorithms can do this
