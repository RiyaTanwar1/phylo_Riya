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

---

## 5. Project Revision: Dataset Change

The original project aimed to analyze whole-genome sequencing (WGS) data from Asian elephants using raw Illumina FASTQ reads. However, due to computational constraints and the scope of this assignment (which focuses on multiple sequence alignment methods), the dataset was revised.

Instead of raw WGS reads, a curated dataset of complete mitochondrial genome sequences was selected. This dataset is more appropriate for evaluating multiple sequence alignment methods and constructing phylogenetic relationships.

---

## 6. Mitochondrial Genome Dataset

A dataset consisting of 10 complete mitochondrial genome sequences (~16–17 kb each) was compiled from NCBI GenBank. Taxa include:

- Elephas maximus  
- Loxodonta africana
- Loxodonta cyclotis  
- Mammuthus primigenius  

All sequences were downloaded in FASTA format and combined into a single multi-FASTA file for alignment.

---

## 7. Alignment Method Selection

The Alignathon study indicates that no single alignment method performs best across all datasets. Performance depends on sequence divergence, dataset size, and structural complexity.

Because this dataset consists of moderately similar mitochondrial genomes and contains 10 sequences, MUSCLE was selected as the alignment method.

---

## 8. Alignment Method: MUSCLE v3.8.31

MUSCLE (Edgar 2004) performs multiple sequence alignment using:

1. K-mer based pairwise distance estimation  
2. Guide tree construction  
3. Progressive alignment  
4. Iterative refinement  

---

### Assumptions

- Sequences are homologous mitochondrial genomes  
- No major structural rearrangements  
- Evolution is primarily substitution-driven  

---

### Limitations

- May not perform as well on highly divergent datasets  
- Structural RNA regions may not be optimally aligned  
- Only one MSA method was tested for this homework  

---

### Command Used

```powershell
& "C:\Users\Riya Tanwar\Downloads\muscle3.8.31_i86win32.exe" -in "elephant_mtDNA_all.fasta" -out "elephant_mtDNA_aligned.fasta"

