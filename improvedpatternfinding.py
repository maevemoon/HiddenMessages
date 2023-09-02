# FOUND A MINIMUM SKEW, BUT K-MER BARELY REPEATS - BUT DnaA CAN BIND TO SLIGHTLY MODIFIED DnaA BOXES AS WELL
# FINDING THE HAMMING DISTANCE (TOTAL NUMBER OF MISMATCHES BETWEEN STRINGS) // FINDING THE APPROXIMATE OCCURENCES OF A PATTERN IN A STRING
# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = [] # initializing list of positions
    #range is modified like this because if is it just length of the bigger text, pattern will keep sliding along with empty letters, adding more to the list of positions
    for i in range(len(Text)-len(Pattern)+1):
        if HammingDistance(Text[i:i+len(Pattern)], Pattern) <= d:
            positions.append(i)
    return positions
# alternative
# def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    p = len(Pattern)
    t = len(Text)
    for i in range(t-p+1):
        if HammingDistance(Text[i:i+p], Pattern) <= d:
            count = count+1
    return count

# Input:  Two strings p and q
# Output: An integer value representing the Hamming Distance between p and q.
def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        # we have to assume that the lengths of the two strings are the same
        # then we check each of their positions and see if its symbol is equal to the other strings
        if len(p) == len(q) and p[i] != q[i]:
            count = count+1
    return count

# EXAMPLE
Text = "ATATGATATGTCTGCTGATG"
Pattern = "ATG"
d = 2

print(ApproximatePatternMatching(Text, Pattern, d))
# RESULT = [0, 2, 5, 7, 9, 11, 14, 17]
# at position 0, Pattern matches text (ATA = ATG with 1 mismatch)
# at position 1, Pattern does not match Text (TAT = ATG with 3 mismatches)
# ...
# at position 7. Pattern matches Text (ATG = ATG)
# at position 9, Pattern matches Text (GTC = ATG with 2 mismatches)

