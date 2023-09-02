# A DIFFERENT APPROACH TO FINDING PATTERNS THAN WEEK 1-2'S FREQUENT WORDS WITH MISMATCHES
# MOTIF = PATTERN THAT APPEARS AT LEAST ONCE IN EACH ONE OF SEVERAL DIFFERENT REGIONS SCATTERED THROUGHOUT THE GENOME
# FIRST, FINDING MOTIFS AND RETURNING A COUNT MATRIX
# Input:  A set of k-mers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count  = {}
    # this is the length of the k-mer, but all are the same length
    # it's being set to the length of the first string in the list of strings
    k = len(Motifs[0])
    # this first loop we are defining the symbols we want as A, C, G, T 
    # a list is then made to set it up to exist
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            # this list has all of these symbols with a value of 0 but now is THE LENGTH of our k-mer
            # basically our list looks like A: 0 0 0 0 0 for 5-mers; 5 zeroes for 5-mers
            count[symbol].append(0)
    t = len(Motifs)
    # now we look for each k-mer in our list of strings, and then for each symbol within said k-mer
    for i in range(t):
        for j in range(k):
            # count[symbol] corresponds to the key of the dictionary count
            # count[symbol][j] corresponds to the position in the list assigned to the key
            # the following line assigns the key (symbol) to a nucleotide in Motifs
            symbol = Motifs[i][j]
            # adds 1 to the position in the list assigned to the key
            count[symbol][j] += 1
    return count

# DIVING INTO PROBABILITY - TRANSFORMING COUNT INTO A PROFILE
# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    # like for count matrix, profile matrix is also a dictionary
    # it is formed the same as the count matrix but then divided by t
    profile = {}
    profile = Count(Motifs)
    k = len(Motifs[0])
    t = len(Motifs)
    # now, for every single element of every single string...
    for i in profile:
        for j in range(k):
            # for our profile matrix with i rows and j columns, it is now redefined as the same but divided by t
            # now each element is divided by t and this makes our profile matrix
            profile[i][j] = profile[i][j]/t
    return profile


# CREATING A CONSENSUS STRING - FINDING THE MOST PROBABLE BASE PER POSITION
# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)
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


# SCORING THE CONSENSUS STRING - THE NUMBER OF UNPOPULAR LETTERS IN THE MOTIF MATRIX
# OUR GOAL IS TO ULTIMATELY FIND THE LOWEST-SCORING CONSENSUS STRING - "MOST PROBABLE"
# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    consensus = Consensus(Motifs)
    score = 0
    for i in range(t):
        for j in range(k):
            # if the i,j position of our Motifs, list of strings, doesn't equal the position of the consensus string, add to score
            if Motifs[i][j] != consensus[j]:
                score += 1
    return score


# GREEDY MOTIF SEARCH - DISASTROUS ALGORITHMS THAT MAY PROVE USEFUL
# FIND THE PROBABILITY THAT Profile(Motifs) GENERATES Text
# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
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

# Input:  A list of k-mers Dna, and integers k and t (where t is the number of k-mers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
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
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j],k,P))
        # once we have our set fully made, we now score it
        # if the score is less than the current score of bestmotifs, we redefine bestmotifs to be the smaller scored motifs
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs