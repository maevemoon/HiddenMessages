# MODIFYING Profile(Motifs) - AN ENTIRE K-MER CAN HAVE A PROBABILITY OF 0 SIMPLY BECAUSE ONE LETTER HAS A PROBABILITY OF 0
# PSEUDOCOUNTS (Laplace's Rule of Succession) ADDS TO SCORING TO REMOVE THIS UNFAIR SCORING
# pseudocounts often amount to adding 1 to each element of Count(Motifs) 

# RETURNING THE PROFILE MATRIX OF MOTIFS WITH PSEUDOCOUNTS
# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ConsensusWithPseudocounts(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    # we are going to build this string one by one - so first define as an empty string
    consensus = ""
    # for every position within every string, let's look for the most frequent symbol
    for j in range(k):
        # we define m as 0 as we're looking at all counts against it (which are >=0)
        # we then define the frequentsymbol as an empty string for now but that will change as we find it
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            # if the current position count we're looking at is greater than m,
            # we make that position our new m, and make our frequentsymbol to that symbol
            # we do this until there is no other symbol in that column with a greater count
            # then we add that frequent symbol to our consensus string, and go on to check all other columns to build the rest
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# finally we find a score using our consensus string and motifs
# our score is calculuated by adding on, column by column, the amount of nucleotides that are NOT present in the consensus string per column
# using this feature we can minimize the score, and thus get the best / most probable consensus string
def Score(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    consensus = ConsensusWithPseudocounts(Motifs)
    score = 0
    for i in range(t):
        for j in range(k):
            # if the i,j position of our Motifs, list of strings, doesn't equal the position of the consensus string, add to score
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score

# now we make a probability function
# using a given string Text and profile matrix we can calculate the probability of a certain k-mer happening
# similarly, we can use this to find the most probable k-mer
def Pr(Text, Profile):
    # we begin by setting a probability value to 1 (certain)
    p = 1
    # now we range through all of text and we keep altering the probability
    # for every character in Text, we take the probability of the nucleotide from our profile matrix at that point
    # we then multiply the probability of that by our current value of p, until we're done
    # eg. Text = AAAAA and A: .1 .2 .3 .4 .5 from profile matrix
    #     p = 1 * .1 * .2 * .3 * .4 * .5 = 0.0012
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

# now let's find a profile-most probable k-mer in a string
# we define a max probability m to be -1 in order for 0 probability strings to still count
# we define a random value x to 0, which will then change to the symbol at the start position of the most probable kmer
def ProfileMostProbableKmer(text, k, profile):
    m = -1
    x = 0
    # as we range through the string we keep redefining our max value until it is the highest probability possible
    # we also keep track of the highest probability's start position in order to then define the kmer
    # the kmer is the window of text starting at x, the most probable start position, going through as long as k is
    for i in range(len(text)-k+1):
        p = Pr(text[i:i+k], profile)
        if p > m:
            m = p
            x = i
    kmer = text[x:x+k]
    return kmer

# GreedyMotifSearch finds the set of motifs across a number of dna sequences that match each other most closely
# we run through the first string and find all kmers
# for all other strings, we find the closest kmers to the original, create a set from it, and then score it
# the lowest score is best - once we find the set that contains them we take the kmers and represent them as a list
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    # we start by setting BestMotifs equal to the first kmer from each DNA string
    # these now count as the best-scoring motifs
    for i in range(0,t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    # now we range through the first dna string for the length of the k-mer
    # we find all the kmers and add them to our list of motifs
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        # this inner loop goes for every set made from the first string, and finds the most probable kmer through all the further strings
        # t = amount of strings
        for j in range(1,t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j],k,P))
        # once we have our set fully made, we now score it
        # if the score is less than the current score of bestmotifs, we redefine bestmotifs to be the smaller scored motifs
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# as we don't want to use probabilities of 0, we use our count function and edit it
# anywhere we would have a count we add 1. therefore nothing would be 0
def CountWithPseudocounts(Motifs):
    pseudocount  = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        pseudocount[symbol] = []
        for j in range(k):
            # instead of count[symbol].append(0) we instead append 1
            # this makes that first initial list start with everything counted to 1
            # then the function proceeds as before
            pseudocount[symbol].append(1)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            pseudocount[symbol][j] += 1
    return pseudocount

# now we're going to do the same thing, but for our profile matrix function
def ProfileWithPseudocounts(Motifs):
    count = CountWithPseudocounts(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    # now, for every single element of every single string...
    for i in range(k):
        total = 0
        for symbol in "ACGT":
            total += count[symbol][i]
        for symbol in "ACGT":
            count[symbol][i] = count[symbol][i]/total
    pseudoprofile = count
    return pseudoprofile

# and now the same with GreedyMotifSearch
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    # we start by setting BestMotifs equal to the first kmer from each DNA string
    # these now count as the best-scoring motifs
    for i in range(0,t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    # now we range through the first dna string for the length of the k-mer
    # we find all the kmers and add them to our list of motifs
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        # this inner loop goes for every set made from the first string, and finds the most probable kmer through all the further strings
        # t = amount of strings
        for j in range(1,t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j],k,P))
        # once we have our set fully made, we now score it
        # if the score is less than the current score of bestmotifs, we redefine bestmotifs to be the smaller scored motifs
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs
# Copy the ten strings occurring in the hyperlinked DosR dataset below.
Dna = ["GCGCCCCGCCCGGACAGCCATGCGCTAACCCTGGCTTCGATGGCGCCGGCTCAGTTAGGGCCGGAAGTCCCCAATGTGGCAGACCTTTCGCCCCTGGCGGACGAATGACCCCAGTGGCCGGGACTTCAGGCCCTATCGGAGGGCTCCGGCGCGGTGGTCGGATTTGTCTGTGGAGGTTACACCCCAATCGCAAGGATGCATTATGACCAGCGAGCTGAGCCTGGTCGCCACTGGAAAGGGGAGCAACATC", "CCGATCGGCATCACTATCGGTCCTGCGGCCGCCCATAGCGCTATATCCGGCTGGTGAAATCAATTGACAACCTTCGACTTTGAGGTGGCCTACGGCGAGGACAAGCCAGGCAAGCCAGCTGCCTCAACGCGCGCCAGTACGGGTCCATCGACCCGCGGCCCACGGGTCAAACGACCCTAGTGTTCGCTACGACGTGGTCGTACCTTCGGCAGCAGATCAGCAATAGCACCCCGACTCGAGGAGGATCCCG", "ACCGTCGATGTGCCCGGTCGCGCCGCGTCCACCTCGGTCATCGACCCCACGATGAGGACGCCATCGGCCGCGACCAAGCCCCGTGAAACTCTGACGGCGTGCTGGCCGGGCTGCGGCACCTGATCACCTTAGGGCACTTGGGCCACCACAACGGGCCGCCGGTCTCGACAGTGGCCACCACCACACAGGTGACTTCCGGCGGGACGTAAGTCCCTAACGCGTCGTTCCGCACGCGGTTAGCTTTGCTGCC", "GGGTCAGGTATATTTATCGCACACTTGGGCACATGACACACAAGCGCCAGAATCCCGGACCGAACCGAGCACCGTGGGTGGGCAGCCTCCATACAGCGATGACCTGATCGATCATCGGCCAGGGCGCCGGGCTTCCAACCGTGGCCGTCTCAGTACCCAGCCTCATTGACCCTTCGACGCATCCACTGCGCGTAAGTCGGCTCAACCCTTTCAAACCGCTGGATTACCGACCGCAGAAAGGGGGCAGGAC", "GTAGGTCAAACCGGGTGTACATACCCGCTCAATCGCCCAGCACTTCGGGCAGATCACCGGGTTTCCCCGGTATCACCAATACTGCCACCAAACACAGCAGGCGGGAAGGGGCGAAAGTCCCTTATCCGACAATAAAACTTCGCTTGTTCGACGCCCGGTTCACCCGATATGCACGGCGCCCAGCCATTCGTGACCGACGTCCCCAGCCCCAAGGCCGAACGACCCTAGGAGCCACGAGCAATTCACAGCG", "CCGCTGGCGACGCTGTTCGCCGGCAGCGTGCGTGACGACTTCGAGCTGCCCGACTACACCTGGTGACCACCGCCGACGGGCACCTCTCCGCCAGGTAGGCACGGTTTGTCGCCGGCAATGTGACCTTTGGGCGCGGTCTTGAGGACCTTCGGCCCCACCCACGAGGCCGCCGCCGGCCGATCGTATGACGTGCAATGTACGCCATAGGGTGCGTGTTACGGCGATTACCTGAAGGCGGCGGTGGTCCGGA", "GGCCAACTGCACCGCGCTCTTGATGACATCGGTGGTCACCATGGTGTCCGGCATGATCAACCTCCGCTGTTCGATATCACCCCGATCTTTCTGAACGGCGGTTGGCAGACAACAGGGTCAATGGTCCCCAAGTGGATCACCGACGGGCGCGGACAAATGGCCCGCGCTTCGGGGACTTCTGTCCCTAGCCCTGGCCACGATGGGCTGGTCGGATCAAAGGCATCCGTTTCCATCGATTAGGAGGCATCAA", "GTACATGTCCAGAGCGAGCCTCAGCTTCTGCGCAGCGACGGAAACTGCCACACTCAAAGCCTACTGGGCGCACGTGTGGCAACGAGTCGATCCACACGAAATGCCGCCGTTGGGCCGCGGACTAGCCGAATTTTCCGGGTGGTGACACAGCCCACATTTGGCATGGGACTTTCGGCCCTGTCCGCGTCCGTGTCGGCCAGACAAGCTTTGGGCATTGGCCACAATCGGGCCACAATCGAAAGCCGAGCAG", "GGCAGCTGTCGGCAACTGTAAGCCATTTCTGGGACTTTGCTGTGAAAAGCTGGGCGATGGTTGTGGACCTGGACGAGCCACCCGTGCGATAGGTGAGATTCATTCTCGCCCTGACGGGTTGCGTCTGTCATCGGTCGATAAGGACTAACGGCCCTCAGGTGGGGACCAACGCCCCTGGGAGATAGCGGTCCCCGCCAGTAACGTACCGCTGAACCGACGGGATGTATCCGCCCCAGCGAAGGAGACGGCG", "TCAGCACCATGACCGCCTGGCCACCAATCGCCCGTAACAAGCGGGACGTCCGCGACGACGCGTGCGCTAGCGCCGTGGCGGTGACAACGACCAGATATGGTCCGAGCACGCGGGCGAACCTCGTGTTCTGGCCTCGGCCAGTTGTGTAGAGCTCATCGCTGTCATCGAGCGATATCCGACCACTGATCCAAGTCGGGGGCTCTGGGGACCGAAGTCCCCGGGCTCGGAGCTATCGGACCTCACGATCACC"]

# EXAMPLE
t = 10
k = 15
Motifs = GreedyMotifSearchWithPseudocounts(Dna, k, t)
print(Motifs)
print(Score(Motifs))

# RESULTS: 
# best set of Motifs
# scoe = 35
# GGACTTCAGGCCCTA
# GGTCAAACGACCCTA
# GGACGTAAGTCCCTA
# GGATTACCGACCGCA
# GGCCGAACGACCCTA
# GGACCTTCGGCCCCA
# GGACTTCTGTCCCTA
# GGACTTTCGGCCCTG
# GGACTAACGGCCCTC
# GGACCGAAGTCCCCG

# TAKING A PROFILE MATRIX Profile CORRESPONDING TO A LIST OF STRINGS Dna AND RETURNS A LIST OF THE Profile-MOST PROBABLE K-MERS IN EACH Dna STRING
# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    motifs = []
    t = len(Dna)
    k = 4
    for i in range(t):
        motif = ProfileMostProbableKmer(Dna[i],k,Profile)
        motifs.append(motif)
    return motifs 

def ProfileMostProbableKmer(text, k, profile):
    m = -1
    x = 0
    # as we range through the string we keep redefining our max value until it is the highest probability possible
    # we also keep track of the highest probability's start position in order to then define the kmer
    # the kmer is the window of text starting at x, the most probable start position, going through as long as k is
    for i in range(len(text)-k+1):
        p = Pr(text[i:i+k], profile)
        if p > m:
            m = p
            x = i
    kmer = text[x:x+k]
    return kmer

def Pr(Text, Profile):
    # we begin by setting a probability value to 1 (certain)
    p = 1
    # now we range through all of text and we keep altering the probability
    # for every character in Text, we take the probability of the nucleotide from our profile matrix at that point
    # we then multiply the probability of that by our current value of p, until we're done
    # eg. Text = AAAAA and A: .1 .2 .3 .4 .5 from profile matrix
    #     p = 1 * .1 * .2 * .3 * .4 * .5 = 0.0012
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p
