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

```



9. Phylogenetic Tree Reconstruction
9.1 Distance-Based Tree (Neighbor-Joining)
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
# Calculate genetic distances using TN93 model
D <- dist.dna(dna, model="TN93")

# Construct Neighbor-Joining tree
tre_nj <- nj(D)

# Improve visualization
tre_nj <- ladderize(tre_nj)

# Plot
plot(tre_nj, cex=0.6)
title("Neighbor-Joining Tree (TN93 distances)")
9.2 Parsimony-Based Tree
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
# Convert alignment for parsimony analysis
dna2 <- as.phyDat(dna)

# Initial tree
tre_ini <- nj(dist.dna(dna, model="raw"))

# Calculate parsimony score
parsimony(tre_ini, dna2)

# Optimize tree under maximum parsimony
tre_pars <- optim.parsimony(tre_ini, dna2)

# Plot
plot(tre_pars, cex=0.6)
title("Maximum Parsimony Tree")
9.3 Comparison of Methods

NJ uses model-corrected distances (TN93)

Parsimony uses character-state changes directly

NJ is generally faster and scalable

Parsimony may be sensitive to long-branch attraction

Both methods assume homologous aligned sequences

Because mitochondrial genomes are relatively small (~16–17 kb) and moderately conserved, both approaches are computationally feasible and appropriate for this dataset.