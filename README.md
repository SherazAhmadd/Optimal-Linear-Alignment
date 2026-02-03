# Optimal-Linear-Alignment
The Myers-Miller algorithm is designed to solve memory and processing bottlenecks in large-scale genomic data. This repository provides a robust **C++ implementation of the Myers-Miller algorithm (1988)** for global sequence alignment. The tool is specifically designed to perform optimal pairwise alignments in **linear space**, addressing the memory constraints inherent in standard dynamic programming approaches when applied to large-scale genomic datasets.

---

## Introduction

In bioinformatics, the comparison of primary biological sequences is fundamental. The traditional **Needleman–Wunsch** algorithm for global alignment requires  
**O(NM)** space complexity to store the scoring matrix, where *N* and *M* represent sequence lengths. For modern genomic sequences exceeding `10^5` base pairs, the memory requirement becomes prohibitive, often surpassing the capacity of standard computing environments.

The **Myers-Miller algorithm** utilizes a divide-and-conquer strategy, based on the principles of the **Hirschberg algorithm**, to achieve optimal alignment using only linear space.

- **Space Complexity:** `O(min(N, M))`, requiring memory proportional only to the sequence length  
- **Time Complexity:** `O(NM)`, maintaining the same quadratic time complexity as standard approaches while reducing the memory footprint by several orders of magnitude

---

## Features

- **Linear Space Efficiency**  
  Optimized for the alignment of long genomic or proteomic sequences without risk of memory overflow.

- **Comprehensive FASTA Support**  
  Integrated parser for biological sequence formats including `.fasta`, `.fna`, `.faa`, and `.fa`.

- **Cross-Platform Normalization**  
  Automatically handles carriage return variations between Windows (`\r\n`) and Unix-based (`\n`) systems to ensure consistent sequence length counts.

- **Quantitative Reporting**  
  Generates detailed alignment output in `.fna` format, including metrics for gap counts, substitutions, and total alignment score.

- **Performance Optimization**  
  Developed in C++ with support for high-level compiler optimizations.

---

## Installation and Setup

### Prerequisites

- A C++ compiler supporting the **g++** standard or later (e.g., `g++ 5.x+`, `clang++`)
- `git` command-line tools

---

### Cloning the Repository

To clone the source code and navigate into the project directory, execute the following commands:

```bash
git clone https://github.com/SherazAhmadd/Optimal-Linear-Alignment.git
cd Optimal-Linear-Alignment
```

### Compilation

For optimal execution speed, it is recommended to compile the source code with the -O3 optimization flag:
```
g++ -O3 myers_miller.cpp -o myers_miller
```
### Usage

The executable requires two sequence files as input arguments. The algorithm will align the first sequence identified within each provided file.
```
./myers_miller test_data/seq1.fna test_data/seq2.fna
```

### Output Data

Upon completion, the program generates an alignment file named:

```
alignment_result.txt
```
This file contains the optimally aligned sequences in a multi-line FASTA format. Summary statistics are also output to the standard console for immediate review.

### Scoring Parameters

The implementation utilizes a discrete global scoring system, which can be modified within the source code constants:

- Match Reward: +2
- Mismatch Penalty: -1
- Gap Penalty: -2

### Contributor / Author
**Name:** Rana Sheraz Ahmad

**Role:** Implementer

**GitHub:** [SherazAhmadd](https://github.com/SherazAhmadd)

## Contact & Issues
Email: ranasheraz.202101902@gcuf.edu.pk 

Issues: [https://github.com/SherazAhmadd/Optimal-Linear-Alignment/issues](https://github.com/SherazAhmadd/Optimal-Linear-Alignment/issues)

## Acknowledgment
Eugene W. Myers, Webb Miller

*[Optimal alignments in linear space](https://doi.org/10.1093/bioinformatics/4.1.11)*

