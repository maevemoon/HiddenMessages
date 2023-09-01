# HiddenMessages

DNA is known to contain repetitive sequences within its strands to which other components can bind, sometimes indicating a notable region. An example of a region like this is the origin of replication in prokaryotic organisms, or oriC, with DnaA boxes as their most abundant repeated sequences.
With information on the frequency and location of repetitive sequences within a genome, scientists can further locate these notable regions. Knowing the locations of such regions is crucial in fields of genetic diagnostics, drug development, genetic engineering, and more.

HiddenMessages is an approximate 1,000 line python project using Motif Discovery, Laplace's Law of Succession, and Gibbs Sampling to identify recognisable DNA sequences within small microbial genomes, such as the oriC of *E. coli*. It also visualises these regions on skew diagrams.

Each subsequent file builds off its predecessor, addressing its challenges and presenting solutions to the gaps identified. Ultimately, each file contains a different approach towards the same goal: finding the most frequently-repeated DNA sequences within a strand of DNA. 

## 1. **patternfinding.py**
The algorithms in this file form the foundation of the HiddenMessages project.
They are introductory, aiming to find and count the most frequently-repeated sequences in DNA strands, or, to put it simply, patterns in strings.

Definitions of the following subroutines:
- PatternCount(Text, Pattern): counts how many times a specified pattern Pattern appears in string Text
- FrequentWords(Text, k): finds and lists the most frequent words of length k that appear in string Text
- FrequencyMap(Text, k): makes a dictionary of the most frequent words of length k that appear in string Text with their count
- ReverseComplement(Pattern): creates a reverse complement of a specified base pattern Pattern
- Reverse(Pattern): reverses a specified base pattern Pattern (e.g. ATG -> GTA)
- Complement(Pattern): complements a specified base pattern Pattern (e.g. ATG -> TAC)
- PatternMatching(Pattern, Genome): lists every starting position for a specified base pattern Pattern within a Dna string Genome
 
## 2. x

## 3. x

## 4. x


This project was supported by the "Bioinformatics for Beginners" course authorised by UC San Diego.
