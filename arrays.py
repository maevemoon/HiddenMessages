# LOCATING THE FORWARD AND REVERSE HALF-STRANDS TO FIND ORI - UTILIZING THE FREQUENCY OF CYTOSINE (due to deanimation) USING SymbolArray, FasterSymbolArray and PatternCount
# a window with fewest occurences of C roughly corresponds to forward half-strand
# a window with most occurences of C roughly corresponds to reverse half-strand

# COUNTING OCCURENCES OF A SYMBOL IN A CIRCULAR GENOME IN EACH WINDOW 
# Input:  Strings Genome and symbol
# Output: SymbolArray(Genome, symbol)
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    # to account for circular genomes, meaning a dna string wraps around itself, extend the genome 
    # since we assume ter is opposite ori and we're looking at one circular strand, the length of our window is n // 2
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(symbol, ExtendedGenome[i:i+(n//2)])
    return array

# A FASTER ALGORITHM: COUNTING OCCURENCES OF A SYMBOL IN A CIRCULAR GENOME BY NOT COUNTING THE SAME CHARACTER TWICE AS THE WINDOW IS MOVED
# Input:  Strings Genome and symbol
# Output: FasterSymbolArray(Genome, symbol)
def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    # we start by using patterncount once
    # this measures how many wanted symbols are in the whole first window
    array[0] = PatternCount(Genome[0:n//2], symbol)
    # this for loop looks at the first and last symbols only
    # we start at 1 and not at 0 because the first thing we do is equal i to previous i's count
    # if we started at 0, we'd be equaling i to -1 i's count which doesn't make sense
    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        # basically, we're not counting again, but just altering the old count
        array[i] = array[i-1]
        # checking if the first symbol of the window before we moved it was what we need
        # if so, we reduce count by 1
        # this removes it from being counted in our current window since it's obviously no longer in our window
        # otherwise we skip this step; count stays the same
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1
        # checking if the symbol that is now at the end of the window is what we need
        # if so, we increase count by 1
        # this is because if it is, it wasn't in the window before so it was not counted
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

def PatternCount(Pattern, Text):
    count = 0
    # "for i in range(n)" iterates for all values of i between 0 and n-1, so we add 1
    # this way it's like we're following the length of the string, not the positions (which start from 0)
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count = count+1
    return count


# KEEPING TRACK OF THE DIFFERENCE BETWEEN THE TOTAL OCCURENCES OF G VS C USING A SKEW ARRAY
# FINDING THE POSITION IN A GENOME WHERE THE SKEW DIAGRAM IS AT A MINIMUM
# Input:  A DNA string Genome
# Output: A list containing all integers i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)
def MinimumSkew(Genome):
    positions = []
    sarray = SkewArray(Genome)
    m = min(sarray)
    for i in range(len(sarray)):
        if sarray[i] == m:
            positions.append(i)
    return positions

# Input:  A String Genome
# Output: SkewArray(Genome)
def SkewArray(Genome):
    skew = [0]
    score = {"A":0, "T":0, "C":-1, "G":1}
    for i in range(len(Genome)):
        # as we traverse the genome we add to our list
        # we add the score of the symbol in the genome position
        # but we also add the skew[i] value so that the score is carried on until it changes by a +-1
        # eg. CAT adding only score(genome[i]) would be -1, 0, 0
        #     but adding skew[i] would keep it -1, -1, -1 (since it's still negative)
        # when appending you're creating a new score which isn't defined yet
        # therefore skew[i] is the previous score you're adding on
        # and you don't have to write skew[i-1]
        skew.append(score[Genome[i]] + skew[i])
    return skew