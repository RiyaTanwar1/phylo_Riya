# Elephant Genomics Project Notebook

Author: Riya Tanwar  
Date Started: February 2026  

---

## 1. Dataset Selection

Dataset: Asian elephant (Elephas maximus) whole-genome sequencing data  
Source: Kappelhof et al., "Genomics Reveals Distinct Evolutionary Lineages in Asian Elephants"

Files selected:

- ERR361996_1.fq.gz
- ERR361996_2.fq.gz

These files represent paired-end Illumina whole-genome sequencing reads from a single Asian elephant individual.

---

## 2. Data Stage

The dataset consists of raw FASTQ sequencing reads.  
Because the reads are unprocessed, quality control is required before downstream analysis.

---

## 3. Quality Control with FastQC

FastQC was run using the following command:

```bash
fastqc ERR361996_1.fq.gz ERR361996_2.fq.gz

---

## 4. Stage of the Dataset

The dataset consists of raw paired-end Illumina sequencing reads for Asian elephants.  

Because the reads are raw, they require quality control before any downstream analysis such as assembly or variant calling.  

**Quality Control:** FastQC was used to check the quality of the raw reads. FastQC evaluates:

- Per-base sequence quality
- GC content
- Adapter contamination
- Sequence duplication levels
- Overrepresented sequences
