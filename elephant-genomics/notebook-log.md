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
```

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

My original project aimed to analyze whole-genome sequencing (WGS) data from Asian elephants using raw Illumina FASTQ reads. However, due to computational constraints and the scope of this assignment (which focuses on multiple sequence alignment methods), the dataset was revised.

Instead of raw WGS reads, I chose to use  complete mitochondrial genome sequences. This dataset is more appropriate for evaluating multiple sequence alignment methods and constructing phylogenetic relationships.

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

- May perform less well on highly divergent datasets  
- Structural RNA regions may not be optimally aligned  
- Only one MSA method was tested for this homework  

---

### Command Used

```powershell
& "C:\Users\Riya Tanwar\Downloads\muscle3.8.31_i86win32.exe" -in "elephant_mtDNA_all.fasta" -out "elephant_mtDNA_aligned.fasta"
```


### 9. Phylogenetic Tree Reconstruction

#### 9.1 Distance-Based Tree (Neighbor-Joining)
Algorithm Used: Neighbor-Joining (NJ)

The Neighbor-Joining algorithm (Saitou & Nei 1987), implemented in the ape package, is a distance-based method that constructs a tree from a matrix of pairwise genetic distances.

Steps:

Compute pairwise genetic distances using dist.dna()

Iteratively join the pair of taxa that minimizes total tree length

Continue clustering until a fully resolved unrooted tree is obtained

In this analysis, distances were calculated using the TN93 substitution model.

Assumptions

Distances accurately reflect evolutionary divergence

Sites evolve independently

Substitution model (TN93) appropriately captures nucleotide evolution

No recombination (reasonable for mitochondrial genomes)

Limitations

Information is reduced to pairwise distances (site-level information lost)

Sensitive to incorrect distance estimation

Produces an unrooted tree

Does not explicitly model character evolution during tree search

Code
#### Calculate genetic distances using TN93 model
D <- dist.dna(dna, model="TN93")

#### Construct Neighbor-Joining tree
tre_nj <- nj(D)

#### Improve visualization
tre_nj <- ladderize(tre_nj)

#### Plot
plot(tre_nj, cex=0.6)
title("Neighbor-Joining Tree (TN93 distances)")

#### 9.2 Parsimony-Based Tree
Algorithm Used: Maximum Parsimony

Maximum Parsimony (implemented in the phangorn package) searches for the tree that minimizes the total number of evolutionary changes required to explain the observed sequence data.

Steps:

Convert DNA alignment to phyDat format

Start with an initial tree (NJ with raw distances)

Optimize tree topology to minimize total character changes

Return tree with smallest parsimony score

Assumptions

Evolution proceeds via the smallest possible number of changes

All substitutions are equally weighted

Homoplasy (convergent evolution) is minimal

Sites evolve independently

Limitations

Sensitive to long-branch attraction

Does not incorporate an explicit substitution model

Can become computationally expensive for large datasets

May be less accurate when substitution rates vary across lineages

Code
#### Convert alignment for parsimony analysis
dna2 <- as.phyDat(dna)

#### Initial tree
tre_ini <- nj(dist.dna(dna, model="raw"))

#### Calculate parsimony score
parsimony(tre_ini, dna2)

#### Optimize tree under maximum parsimony
tre_pars <- optim.parsimony(tre_ini, dna2)

#### Plot
plot(tre_pars, cex=0.6)
title("Maximum Parsimony Tree")

#### 9.3 Comparison of Methods

NJ uses model-corrected distances (TN93)

Parsimony uses character-state changes directly

NJ is generally faster and scalable

Parsimony may be sensitive to long-branch attraction

Both methods assume homologous aligned sequences

Because mitochondrial genomes are relatively small (~16–17 kb) and moderately conserved, both approaches are computationally feasible and appropriate for this dataset.

### 10. Maximum Likelihood Phylogenetic Analysis using IQ-TREE

Description:
IQ-TREE uses a maximum likelihood approach to infer phylogenetic trees. It searches for the tree topology that maximizes the probability of the observed sequence data under a model of sequence evolution. It uses efficient tree search algorithms and ModelFinder to select the best substitution model.

Assumptions: Sequences are homologous and correctly aligned Sites evolve independently. Evolution follows a substitution model (e.g., GTR). The process is stationary, reversible, and homogeneous. 

Limitations: Sensitive to alignment errors, computationally intensive for large datasets, model misspecification can bias results, may get stuck in local optima.

#### Run IQ-TREE (version 3)

$iqtree = ".\iqtree-3.0.1-Windows\iqtree-3.0.1-Windows\bin\iqtree3.exe"

& $iqtree `
  -s .\data\elephant_mtDNA_aligned.fasta `
  -m MFP `
  -bb 1000 `
  -alrt 1000


### 11. Bayesian Phylogenetic Analysis using MrBayes

Description:
Bayesian phylogenetic inference estimates the posterior probability distribution of trees given the observed sequence data and a substitution model. MrBayes implements a Markov Chain Monte Carlo (MCMC) approach to sample trees according to their posterior probability. This method provides statistical support for clades in the form of posterior probabilities.

Assumptions:

Sequences are homologous and correctly aligned
Sites evolve independently
Evolutionary process is stationary, reversible, and homogeneous
The substitution model (here HKY + gamma) appropriately describes nucleotide evolution

Limitations:

Computationally intensive for large datasets
Sensitive to alignment errors and model misspecification
Requires careful assessment of MCMC convergence and sufficient burn-in
MrBayes Command Block

begin mrbayes;
    set autoclose=yes;
    prset brlenspr=unconstrained:exp(10.0);
    prset shapepr=exp(1.0);
    prset tratiopr=beta(1.0,1.0);
    prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
    lset nst=2 rates=gamma ngammacat=4;
    mcmcp ngen=1000000 samplefreq=50 printfreq=500 nruns=2 nchains=4 savebrlens=yes;
    outgroup Anacystis_nidulans;
    mcmc;
    sumt burnin=12500;
end;

Explanation of Parameters
ngen = 1,000,000: Increases the number of MCMC generations for better convergence of posterior distributions.
samplefreq = 50: Samples every 50 generations to reduce output file size while retaining sufficient resolution.
printfreq = 500: Prints MCMC progress less frequently to avoid clutter.
nruns = 2: Two independent MCMC runs allow comparison of results to assess convergence.
nchains = 4: Uses one cold chain and three heated chains to improve mixing and escape local optima.
sumt burnin = 12,500: Discards the first 25% of sampled trees as burn-in to ensure only post-convergence trees are summarized.
outgroup = Anacystis_nidulans: Rooting the tree with an outgroup allows proper directionality in the phylogeny.

Output:

A majority-rule consensus tree with posterior probabilities for each clade.
Estimated branch lengths and support values for interpretation of evolutionary relationships among Asian elephant mitochondrial genomes.

Rationale:
Bayesian inference complements the NJ and parsimony analyses, providing a probabilistic framework that accounts for uncertainty in tree topology and branch lengths. Given the small size and moderate divergence of the mitochondrial genome dataset, this approach is computationally feasible and yields robust posterior support for clades.


### 12. ASTRAL Analysis

ASTRAL (Accurate Species TRee ALgorithm) is a summary coalescent method that estimates a species tree from a set of gene trees. Instead of concatenating sequences, ASTRAL uses quartet-based inference, this means it evaluates relationships among all sets of four taxa and finds the species tree that agrees with the largest number of gene tree quartets.

This approach is grounded in coalescent theory, which models how gene lineages trace back through time and accounts for incomplete lineage sorting (ILS)—a common source of discordance between gene trees and species trees.

Assumptions
Gene trees are correctly estimated (or at least not heavily biased)
Discordance among gene trees is primarily due to incomplete lineage sorting (ILS)
No significant effects from:
Horizontal gene transfer
Hybridization
Gene duplication/loss (ASTRAL assumes one copy per species per gene)
Loci are independent

Strengths
Handles large genomic datasets efficiently
Statistically consistent under the multispecies coalescent model
More robust than concatenation when ILS is high
Does not require sequence alignments, only gene trees

Limitations
Sensitive to gene tree estimation error
Cannot explicitly model hybridization or introgression
Assumes no paralogs (single-copy genes)
Support values are local posterior probabilities, not traditional bootstrap

Commands Used 
Step 1: Download ASTRAL
git clone https://github.com/smirarab/ASTRAL.git
cd ASTRAL

Step 2: Run ASTRAL on toy dataset
java -jar astral.5.7.8.jar \
 -i song_mammals.424.gene.tre \
 -o song-astral.tre

Step 3: Reading the Tree in R
library(ape)

tre <- read.tree("song-astral.tre")
plot(tre)

nodelabels(text = tre$node.label)