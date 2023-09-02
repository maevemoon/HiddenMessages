# Hidden Messages

DNA is known to contain repetitive sequences within its strands to which other components can bind, sometimes indicating a notable region. An example of a region like this is the origin of replication in prokaryotic organisms, or oriC, with DnaA boxes as their most abundant repeated sequences.
With information on the frequency and location of repetitive sequences within a genome, scientists can further locate these notable regions. Knowing the locations of such regions is crucial in fields of genetic diagnostics, drug development, genetic engineering, and more.

Hidden Messages is an approximate 1,000 line python project using Motif Discovery, Laplace's Law of Succession, and Gibbs Sampling to identify recognisable DNA sequences and features within small microbial genomes, such as the oriC of *E. coli*. 

Each subsequent file builds off its predecessor, addressing its challenges and presenting solutions to the gaps identified. Ultimately, each file contains a different approach towards the same goal: finding the most frequently-repeated DNA sequences within a strand of DNA. 

## 1. patternfinding.py
The algorithms in this file form the foundation of the Hidden Messages project.
They are introductory, aiming to find and count the most frequently-repeated sequences in DNA strands, or, to put it simply, patterns in strings.

Definitions of the following subroutines:
- PatternCount(Text, Pattern): counts how many times a specified pattern Pattern appears in string Text
- FrequentWords(Text, k): finds and lists the most frequent words of length k that appear in string Text
- FrequencyMap(Text, k): makes a dictionary of the most frequent words of length k that appear in string Text with their count
- ReverseComplement(Pattern): creates a reverse complement of a specified base pattern Pattern
- Reverse(Pattern): reverses a specified base pattern Pattern (e.g. ATG -> GTA)
- Complement(Pattern): complements a specified base pattern Pattern (e.g. ATG -> TAC)
- PatternMatching(Pattern, Genome): lists every starting position for a specified base pattern Pattern within a Dna string Genome
 
## 2. arrays.py
The algorithms in this file take on a different approach at localising oriC. Rather than searching for frequently repeated sequences, they help locate the forward and reverse half-strands by calculating the G/C content in a genome, the composition of which changes across the genome as a result of cytosine undergoing deanimation. A skew diagram is also used, highlighting the minimum skew and therefore indicating the location of oriC. This file also takes into account the circular genomes of prokaryotes, as the file before it assumed genomes to be linear.

Definitions of the following subroutines:
- SymbolArray(Genome, symbol): counts the occurences of a symbol symbol, using PatternCount, in a circular genome Genome in each moving window of length n (n = len(Genome))
- FasterSymbolArray(Genome, symbol): a faster, more effective algorithm than SymbolArray; the same character is not counted more than once if it appears in a different window
- SkewArray(Genome): cumultatively scores a genome based on G/C content across the length of the genome
- MinimumSkew(Genome): lists all positions within the genome where SkewArray reaches a minimum
  
## 3. improvedpatternfinding.py
These two algorithms directly build off patternfinding.py by allowing for mismatches between strings, as a common occurence in DNA. The amount of these mismatches can be referred to as the hamming distance; the number of positions in which the two bases are different when comparing two strings of DNA. DnaA can bind to slightly modified DnaA boxes in the oriC, so these algorithms take this into account and return the starting positions of all occurences of a known Pattern with at most d mismatches.

Definitions of the following subroutines:
- ApproximatePatternMatching(Text, Pattern, d): lists every starting position for a specified pattern Pattern within a string Text with at most d mismatches
- HammingDistance(p, q): returns an integer value representing the hamming distance between two strings p and q

## 4. motifs.py
The algorithms in this file provide a completely different approach into finding and counting the most frequently-repeated sequences in DNA strands. Unlike the two PatternMatching algorithms, it is assumed the pattern of the most frequent k-mer in DNA is not known. Unlike FrequentWords, it uses Motif Discovery to gradually create the most probable consensus string, acknowledging mismatches, to act as the most frequent k-mer in a DNA string.

Definitions of the following subroutines:
- Count(Motifs): given a collection of motifs Motifs, counts the occurrences of each nucleotide in Motifs at each position i in a string of Motifs of length j. Returns this count as a dictionary
- Profile(Motifs): given a collection of motifs Motifs, computes the count, and then finds the fractional probability of encountering a given nucleotide in Motifs at each position i in a string of Motifs of length j. Returns this profile as a dictionary
- Consensus(Motifs): given a collection of motifs Motifs, finds the most probable sequence of nucleotides in Motifs for each given position. Returns this consensus sequence as a string
- Score(Motifs): computes the consensus string consensus, and then summation counts all mismatches of this consensus string to the strings in Motifs. Returns this count as an integer
- Pr(Text, Profile): given string Text and the probabilities in Profile(Motifs), finds the total probability (independent events) of encountering that string. Returns this probability as a float
- ProfileMostProbableKmer(Text, k, Profile) : given string Text, the length of the k-mer k and probabilities in Profile(Motifs), generate the total probabilities of finding a k-mer in Text and find the k-mer with the highest probability. Return this k-mer as a string

## 5. pseudocounts.py
The algorithms in this file are improved from motifs.py by using Laplace's Rule of Succession to remove unfair scoring to the probabilities calculated previously. An entire k-mer can have a probability of 0 simply because one base has a probability of 0. Pseudocounts are added into these algorithms in order to build a profile where all probabilities are bigger than 0. 

Definitions of the following subroutines:
- GreedyMotifSearch(Dna, k, t): finds the set of motifs Motifs across t strings that match each other most closely, scores it, and returns a list of the k-mers with the lowest scores
- CountWithPseudocounts(Motifs): given a collection of motifs Motifs, counts the occurrences of each nucleotide in Motifs at each position i in a string of Motifs of length j while adding pseudocounts. Returns this count as a dictionary
- ProfileWithPseudocounts(Motifs): given a collection of motifs Motifs, computes the count while adding pseudocounts, and then finds the fractional probability of encountering a given nucleotide in Motifs at each position i in a string of Motifs of length j. Returns this profile as a dictionary
- ConsensusWithPseudocounts(Motifs): given a collection of motifs Motifs, finds the most probable sequence of nucleotides in Motifs for each given position while adding pseudocounts. Returns this consensus sequence as a string
- GreedyMotifSearchWithPseudocounts(Dna, k, t): finds the set of motifs across t strings that match each other most closely while adding pseudocounts, scores it, and returns a list of the k-mers with the lowest scores

## 6. gibbssampling.py
The final few algorithms for the Hidden Messages project using Gibbs Sampling, ultimately combining everything written previously to fairly find, count, and consider the most frequently-repeated sequences in DNA strands. It begins as a function RandomMotifs(Dna, k, t), before it is improved one final time into GibbsSampler(Dna, k, t, n) which can overcome hurdles around the local optima and potentially reach the global optimum.

Definitions of the following subroutines:
- RandomMotifs(Dna, k, t): finds randomly-chosen motifs Motifs at each run
- GibbsSampler(Dna, k, t, n): begin by finding randomly-chosen motifs Motifs, and with every subsequent run, randomly change one string from a a list of t sitrings, change one motif to another by rolling a weighted die based on their probabilities. Return the set of Motifs and the total score
- Normalize(Probabilities): rescales a collection of probabilities so they sum to 1
- WeightedDie(Probabilities): simulates rolling a die so that the probability of rolling the i-th side of the die corresponds to the probability of the i-th k-mer in the list
- ProfileGeneratedString(Text, profile, k): after normalizing and rolling a die, randomly chooses a k-mer from a string Text based on a profile matrix Profile


This project was supported by the "Bioinformatics for Beginners" course authorised by UC San Diego.
