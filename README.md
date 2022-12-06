# GenomicToolBox
Library for various tasks in Bioinformatics. Specifically, read simulation, GO enrichment analysis, detection of exon skipping events and classification of reads in bam files.

For each provided jar file, print the exact usage via jar -jar <jarfile> --help.


### Exon Skipping Events

One form of alternative splicing is exon skipping. An exon-skipping splicing event is a tuple (gene, intron-start, intron-end) and is defined by (at least) two transcripts: wildtype (WT) and spliced variant (SV) of the same gene, an intron in SV with start and end corresponding to an exon end and exon start in WT, respectively, and the SV-intron spans at least one exon in WT. For any exon-skip event there may be several WT-s and several SV-s, and there may be several sets of skipped exons (see figure below; this is one exon-skipping event).

![Alt text](images/exon_skipping_event.png?raw=true "Title")

The program extracts exon-skipping events between coding sequences (CDS) for the provided gtf file and writes the results into a tsv file.


### Simple paired-end read simulator

The results of RNA-seq experiments are the sequenced reads in FASTQ format (see: https://en.wikipedia.org/wiki/FASTQ_format).

For paired-end experiments the experiment yields in two files containing the read pairs in the same order. The simulator writes three files, two for the simulated paired-end sequences in FASTQ format (fw.fastq, rw.fastq, one FASTQ file for the first fragment, one for the second, and a tab separated file (read.mappinginfo) containing details about each simulated read.


### Bam Feature Extraction

Extracts information on all mapped read pairs. 

Provides output for all genic or intergenic reads pairs.

Genic are all transcriptomic, merged-transcriptomic and intronic read pairs (see below), i.e. all read pairs where the read pair is contained within the gene (all mapped bases of the reads are ≥ gene start and < gene end for a non-end exlusive gene representation).

Intergenic are all other read pairs except they contain at least one gene completely, i.e. first mapped base of read ≤ start of gene and last mapped base of read >= of gene. The read pair features annotated are:

- Mismatch count
- Clipping size
- Split count
- Gene count: The number of genes that get annotated by the program for this read pair (i.e. the number of genes with the highest match-plausibility level).

- Gene distance: Distance to the next sense gene (only for intergenic read pairs).
- Antisense: would it match a gene if the read pair was on the other strand (true/false; only if no gene is matched).
- pcr-index (format: pcrindex: ${pcrindex}) the number of read-pairs mapped exactly to the same genomic region as this read pair so far.

A read pair does not match any gene if reads are mapped to different chromosomes. A read pair is split inconsistent (annotate split-inconsistent: true) if there is at least one read base in one of the reads that is within a split of the other read (skip all other annotations in these cases).

A gene is matched, if the read pair is fully contained in the gene body disregarding introns. A read pair matches a transcript, if it is contained in the transcript and intron boundaries are consistent. It matches the merged transcript, if each read base is part of some exon.

A read pair is intergenic if no gene is contained between its first base and its last base.


### GO enrichment analysis

Analyses the Gene Ontology (GO) DAG overlap properties and performs enrichment analysis on the provided data in obo file format.
