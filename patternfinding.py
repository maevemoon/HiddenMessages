# SLIDING A PATTERN THROUGH A DNA STRAND - HOW MANY TIMES DOES THE STRAND CONTAIN THIS PATTERN?
# the first window starts at the first character, but the final window starts at position n-k (len(Text) - len(Pattern))
def PatternCount(Text, Pattern):
    count = 0
    # make sure to add 1 as the last element in a range isn't included as part of the range, but rather marks where it ends
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count 


# FINDING MOST FREQUENT K-MERS IN A DNA STRING & MAKING A FREQUENCY MAP, COUNTING OCCURENCES OF EACH POSSIBILITY 
# Input:  A string Text and an integer k
# Output: A list containing all most frequent k-mers in Text
def FrequentWords(Text, k):
    # in order to put our most frequent k-mers somewhere, we have to make a list
    # we redefine 'freq' to be our frequency map
    # then we define 'm' to be the maximum value in that frequency map
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    # time to search for the new variable 'key' in the program
    for key in freq:
        # if 'key' has the max number of patterns, define 'pattern' with it
        if freq[key] == m:
            pattern = key
            # this adds the pattern to our dictionary
            words.append(pattern)
    # time to return our new dictionary, now with all max values in the genome
    return words

def FrequencyMap(Text, k):
    # generating an empty "dictionary" of our patterns
    freq = {}
    # defining n before we start the for loop for ease to not have to retype
    n = len(Text)
    # using PatternCount
    for i in range(n-k+1):
        # this first part of the loop "scans" the genome and locates new patterns
        # a pattern of length k is created for each new window in the genome
        # this pattern is then assigned a frequency of 0 (we just found it, not counted it)
        # this is good, because we don't have to write out the patterns ourselves
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
        for i in range(n-k+1):
            # this second for loop scans the genome again
            # count/frequency increases every time one of our patterns is found
            if Text[i:i+k] == Pattern:
                freq[Pattern] = freq[Pattern]+1
    return freq

# EXAMPLE: set Text equal to the Vibrio cholerae oriC and k equal to 10, and print the result of calling FrequentWords on Text and k
Text = "ATCAATGATCAACGTAAGCTTCTAAGCATGATCAAGGTGCTCACACAGTTTATCCACAACCTGAGTGGATGACATCAAGATAGGTCGTTGTATCTCCTTCCTCTCGTACTCTCATGACCACGGAAAGATGATCAAGAGAGGATGATTTCTTGGCCATATCGCAATGAATACTTGTGACTTGTGCTTCCAATTGACATCTTCAGCGCCATATTGCGCTGGCCAAGGTGACGGAGCGGGATTACGAAAGCATGATCATGGCTGTTGTTCTGTTTATCTTGTTTTGACTGAGACTTGTTAGGATAGACGGTTTTTCATCACTGACTAGCCAAAGCCTTACTCTGCCTGACATCGACCGTAAATTGATAATGAATTTACATGCTTCCGCGACGATTTACCTCTTGATCATCGATCCGATTGAAGATCTTCAATTGTTAATTCTCTTGCCTCGACTCATAGCCATGATGAGCTCTTGATCATGTTTCCTTAACCCTCTATTTTTTACGGAAGAATGATCAAGCTGCTGCTCTTGATCATCGTTTC"
k = 10
print(FrequentWords(Text, k))



# DNA CONTAINS TWO STRANDS: THE TRANSCRIBED STRAND, AND THE COMPLEMENTARY STRAND - WHAT ARE THE MOST FREQUENT K-MERS TAKING INTO ACCOUNT BOTH STRANDS?

# FIRST, MAKING A REVERSE COMPLEMENT OF DNA
# Input:  A DNA string Pattern
# Output: The reverse complement of Pattern
def ReverseComplement(Pattern):   
    # when we pass a variable to a function, a copy of it is made
    # so it's fine if we define a variable twice to a function; it won't change
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern

def Reverse(Pattern):
    # we start with an empty string to build the strand
    rev = ""
    # in this for loop we are looking for all of the characters in a pattern
    # we go through each character, adding it on to our empty string
    for char in Pattern:
        # order matters when adding strings:
        # since we're reversing, we add a new character before the previous
        rev = char + rev
    return rev

def Complement(Pattern):
    # once again, we start with an empty string so we can build the strand
    comp = ""
    # this command replaces the letter from our pattern to another; the complement
    basepairs = {"A":"T", "T":"A", "C":"G", "G":"C"}
    for char in Pattern:
        comp = comp + basepairs.get(char)
    return comp



# WHILE PatternCount RETURNS A COUNT OF HOW MANY TIMES A PATTERN WAS FOUND IN A DNA STRAND, PatternMatching RETURNS THE STARTING POSITIONS OF THE PATTERN IN A GENOME
# FIND ALL START POSITIONS OF A PATTERN IN A DNA STRAND
# perhaps multiple occurences of a pattern in a string of dna repeat throughout the entire genome, not just the ori region
def PatternMatching(Pattern, Genome):
    positions = []
    g = len(Genome)
    p = len(Pattern)
    for i in range(g-p+1):
        if Genome[i:i+p] == Pattern:
            # this will add the starting position, i, to our list
            positions.append(i)
    return positions

# EXAMPLE: call PatternMatching with Pattern equal to "CTTGATCAT" and Genome equal to v_cholerae, and store the output as a variable called positions
# Pattern = "CTTGATCAT"
# Genome = // IMPORT HERE // 
# positions = PatternMatching(Pattern, Genome)
# print(positions)