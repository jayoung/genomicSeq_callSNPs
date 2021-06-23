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

