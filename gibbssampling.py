# GIBBS SAMPLING

# IMPLEMENTING A RANDOM NUMBER GENERATOR - FINDING RANDOMLY-CHOSEN K-MERS TO GENERATE NEW COLLECTIONS OF K-MERS
import random
# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
def RandomMotifs(Dna, k, t):
    randomm = []
    t = len(Dna)
    l = len(Dna[0])
    # time to go through all strings of Dna and pick a random k-mer per string
    for i in range(t):
        # we are picking a random number for the length of the first string of Dna, minus the length of the k-mer
        # this way we don't choose the first position of the random k-mer to be, for example, the last symbol of the string
        r = random.randint(1, l-k)
        randomm.append(Dna[i][r:r+k])
    return randomm

# RandomizedMotifSearch(Dna, k, t) MAY CHANGE ALL K-MERS IN ONE STEP, ULTIMATELY MOVING FROM HIGHER TO LOWER-SCORING MOTIFS
# THIS MEANS SOME CORRECT MOTIFS MAY BE DISCARDED 
# GibbsSampler(Dna, k, t, n) IS AN IMPROVED, CAUTIOUS ALGORITHM WHICH DISCARDS A SINGLE K-MER ONLY AT EACH STEP , ULTIMATELY MOVING FROM LOWER TO HIGHER-SCORING MOTIFS
# first, randomly choose a string from a list of t strings
# then, randomly choose the most probable k-mer from that string after normalizing the probabilities of every k-mer within it, and rolling a weighted die

# but why not choose the k-mer with the highest probability immediately? why roll dies?
# by allowing to take worse-than-optimal choices every now and then, but placing preference on better choices, you give yourself the opportunity to overcome hills around local optima to potentially reach the global optimum
def GibbsSampler(Dna, k, t, n):
    BestMotifs = [] 
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for j in range(1, n):
        # first choice: random string from list of t strings
        i = random.randint(0,t-1)
        ReducedMotifs = []
        for j in range(0,t):
            if j != i:
                ReducedMotifs.append(Motifs[j])
        Profile = ProfileWithPseudocounts(ReducedMotifs)
        # profile-generated k-mer in the i-th string
        Motif_i = ProfileGeneratedString(Dna[i], Profile, k)
        Motifs[i] = Motif_i
        if Score(Motifs) < Score(BestMotifs):
                BestMotifs = Motifs
    return BestMotifs

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    profile = {} 
    c = CountWithPseudocounts(Motifs)
    for n in 'ACGT':
        p = []
        for i in range(0,k):
            p.append(c[n][i]/(t+4))
        profile[n] = p
    return profile

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = {} 
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
             count[symbol].append(1)
    for i in range(t):
        for j in range(k):
             symbol = Motifs[i][j]
             count[symbol][j] += 1
    return count 

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    # insert your code here
    p = 1
    for i in range(0,len(Text)):
        p *= Profile[Text[i]][i]
    return p

# RESCALING A COLLECTION OF PROBABILITIES SO THAT THEY SUM TO 1
# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    result = {}
    sum = 0
    for m in Probabilities:
        sum += Probabilities[m]
    for n in Probabilities:
        result[n]= Probabilities[n]/sum
    return result  

# SIMULATING ROLLING A DIE SO THAT THE PROBABILITY OF ROLLING THE i-TH SIDE OF THE DIE CORRESPONDS TO THE PROBABILITY OF THE i-TH K-MER IN THE LIST
# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    count = 0
    p = random.uniform(0,1)
    for keys,values in Probabilities.items():
        count = count+values
        if p < count:
        
            return keys
        
# AFTER STIMULATING A WEIGHTED DIE ROLL OVER A COLLECTION OF PROBABILITIES OF STRINGS
# RANDOMLY CHOOSE A K-MER FROM A STRING Text BASED ON A PROFILE MATRIX profile
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    # your code here
    n = len(Text)
    probabilities = {} 
    # ranger over all possible k-mers in Text, computing the probability of each one and placing this probability into a dictionary
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
        # normalize these probabilities and return the result of rolling a weighted die over this dictionary to product a k-mer
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers
def Score(Motifs):
    # Insert code here
    k = len(Motifs[0])
    t = len(Motifs)
    cs = ConsensusWithPseudocounts(Motifs)
    score = 0
    for j in range(0,k):
        for i in range(0,t):
            if Motifs[i][j] != cs[j]:
                score += 1
    return score

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs
def ConsensusWithPseudocounts(Motifs):
    k = len(Motifs[0])
    count = CountWithPseudocounts(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus    

# EXAMPLE
k = 8
t = 5
N = 100
Dna = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA","GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG","TAGTACCGAGACCGAAAGAAGTATACAGGCGT","TAGATCAAGTTTCAGGTGCACGTCGGTGAACC","AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

BestMotifs = GibbsSampler(Dna, k, t, N)
print(BestMotifs)
print(Score(BestMotifs))

# RESULTS
# run 1
# ['GTGTTCAG', 'GGGCGAGG', 'CCGAAAGA', 'GTGCACGT', 'GTGCAATG']
# 13
# run 2
# ['GTTCAGTA', 'GTGTAAGT', 'ATACAGGC', 'GGTGAACC', 'AATCCACC']
# 16
# run 3
# ['TGTTCAGT', 'TGCCAAGG', 'TATACAGG', 'ACGTCGGT', 'CGTGCAAT']
# 14
# run 4
# ['GTTCAGTA', 'GGTATGTG', 'CGAGACCG', 'TTTCAGGT', 'GTTGGCCT']
# 18 
# ...
