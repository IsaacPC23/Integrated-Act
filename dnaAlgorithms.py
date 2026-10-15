# Algorithms and data structures used to analyze the SARS-CoV-2 genome.
# Contains: FASTA reading, KMP string matching, palindrome search, longest common substring and the codon tables.
# Dictionary that maps every amino acid to all of its possible DNA codons.
# A dictionary is used because it answers which codons encode the amino acid in O(1) average time.
AMINO_TO_CODONS = {
    "A": ["GCT", "GCC", "GCA", "GCG"],
    "R": ["CGT", "CGC", "CGA", "CGG", "AGA", "AGG"],
    "N": ["AAT", "AAC"],
    "D": ["GAT", "GAC"],
    "C": ["TGT", "TGC"],
    "Q": ["CAA", "CAG"],
    "E": ["GAA", "GAG"],
    "G": ["GGT", "GGC", "GGA", "GGG"],
    "H": ["CAT", "CAC"],
    "I": ["ATT", "ATC", "ATA"],
    "M": ["ATG"],
    "L": ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"],
    "K": ["AAA", "AAG"],
    "F": ["TTT", "TTC"],
    "P": ["CCT", "CCC", "CCA", "CCG"],
    "S": ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"],
    "T": ["ACT", "ACC", "ACA", "ACG"],
    "W": ["TGG"],
    "Y": ["TAT", "TAC"],
    "V": ["GTT", "GTC", "GTA", "GTG"],
    "*": ["TAA", "TGA", "TAG"]
}

# buildCodonToAmino: builds the inverse table of AMINO_TO_CODONS, so one codon can be translated into its amino acid in O(1) average time.
# Returns: dictionary codon -> amino acid letter (* means stop codon).
# Complexity: O(64), constant, it only walks the 64 codons of the table.
def buildCodonToAmino():
    table = {}
    for amino in AMINO_TO_CODONS:
        for codon in AMINO_TO_CODONS[amino]:
            table[codon] = amino
    return table

CODON_TO_AMINO = buildCodonToAmino()

# readFasta: reads a file in FASTA format.
# filePath: path of the file to read.
# Returns: a tuple (header, sequence); the first line is the description and
# the remaining lines are joined into a single continuous string.
# Complexity: O(n), where n is the number of characters in the file.
def readFasta(filePath):
    header = ""
    pieces = []
    inputFile = open(filePath, "r")
    lines = inputFile.readlines()
    inputFile.close()
    if len(lines) > 0:
        header = lines[0].strip()
    for index in range(1, len(lines)):
        pieces.append(lines[index].strip())
    sequence = "".join(pieces)
    return header, sequence

# readFastaMultiple: reads a FASTA file that contains several records, like the protein file, where every record starts with its own ">" line.
# Returns: list of tuples in the order they appear.
# Complexity: O(n), where n is the number of characters in the file.
def readFastaMultiple(filePath):
    records = []
    currentName = ""
    currentPieces = []
    inputFile = open(filePath, "r")
    lines = inputFile.readlines()
    inputFile.close()
    for line in lines:
        cleanLine = line.strip()
        if len(cleanLine) > 0 and cleanLine[0] == ">":
            if len(currentName) > 0:
                records.append((currentName, "".join(currentPieces)))
            currentName = cleanLine[1:]
            currentPieces = []
        else:
            currentPieces.append(cleanLine)
    if len(currentName) > 0:
        records.append((currentName, "".join(currentPieces)))
    return records

# buildPrefixTable: builds the prefix/failure table used by the KMP algorithm.
# pattern: string whose borders are computed.
# Returns: list where position i holds the length of the longest proper prefix
# of pattern[0..i] that is also a suffix of pattern[0..i].
# complexity: O(m), m is the length of the pattern.
def buildPrefixTable(pattern):
    table = [0] * len(pattern)
    length = 0
    index = 1
    while index < len(pattern):
        if pattern[index] == pattern[length]:
            length = length + 1
            table[index] = length
            index = index + 1
        elif length > 0:
            length = table[length - 1]
        else:
            table[index] = 0
            index = index + 1
    return table

# kmpSearch: searches every occurrence of a pattern inside a text with the KMP algorithm.
# pattern: string that is searched
# Returns: list with the starting indexes (0 based) of all the occurrences.
# Complexity: O(n+m), where n is the text length and m the pattern length.
def kmpSearch(text, pattern):
    positions = []
    if len(pattern) == 0 or len(text) < len(pattern):
        return positions
    table = buildPrefixTable(pattern)
    textIndex = 0
    patternIndex = 0
    while textIndex < len(text):
        if text[textIndex] == pattern[patternIndex]:
            textIndex = textIndex + 1
            patternIndex = patternIndex + 1
            if patternIndex == len(pattern):
                positions.append(textIndex - patternIndex)
                patternIndex = table[patternIndex - 1]
        elif patternIndex > 0:
            patternIndex = table[patternIndex - 1]
        else:
            textIndex = textIndex + 1
    return positions

# longestPalindrome: finds the longest palindromic substring expanding from every possible center of the sequence. Every position is tried as the 
# center of a palindrome of odd length, and every space between two positions as the center of a palindrome of even length.
# sequence: string to analyze.
# Returns: tuple (start, length) with the position and size of the longest
# palindrome; (0, 0) when the sequence is empty.
# Complexity: O(n*k) where k is the size of the longest palindrome found /  O(n * n) in the worst case, for example when every letter is the same.
def longestPalindrome(sequence):
    bestStart = 0
    bestLength = 0
    for center in range(len(sequence)):
        left = center
        right = center
        while left >= 0 and right < len(sequence) and sequence[left] == sequence[right]:
            left = left - 1
            right = right + 1
        if (right - left - 1) > bestLength:
            bestLength = right - left - 1
            bestStart = left + 1
        left = center
        right = center + 1
        while left >= 0 and right < len(sequence) and sequence[left] == sequence[right]:
            left = left - 1
            right = right + 1
        if (right - left - 1) > bestLength:
            bestLength = right - left - 1
            bestStart = left + 1
    return bestStart, bestLength

# longestCommonSubstring: finds the longest substring shared by two sequences using dynamic programming with tabulation. Only two
# rows of the table are kept, because row i only depends on row i - 1.
# Returns: tuple (startInFirst, startInSecond, length) of the longest common
# substring; lengths of 0 mean that the sequences share nothing.
# Complexity: O(n*m) time and O(m) extra memory, with n and m the lengths.
def longestCommonSubstring(first, second):
    previousRow = [0] * (len(second) + 1)
    bestLength = 0
    bestEndInFirst = 0
    bestEndInSecond = 0
    for i in range(1, len(first) + 1):
        currentRow = [0] * (len(second) + 1)
        for j in range(1, len(second) + 1):
            if first[i - 1] == second[j - 1]:
                currentRow[j] = previousRow[j - 1] + 1
                if currentRow[j] > bestLength:
                    bestLength = currentRow[j]
                    bestEndInFirst = i
                    bestEndInSecond = j
        previousRow = currentRow
    return bestEndInFirst - bestLength, bestEndInSecond - bestLength, bestLength

# translateSequence: translates a DNA sequence into its amino acid sequence, reading it codon by codon from a given reading frame.
# sequence: DNA sequence formed by the letters A, C, G and T.
# frame: 0, 1 or 2, the position where the first codon starts.
# Returns: string of amino acid letters, "*" for stop codons and "?" for codons that are not in the table.
# Complexity: O(n), where n is the length of the sequence.
def translateSequence(sequence, frame):
    aminoAcids = []
    for index in range(frame, len(sequence) - 2, 3):
        codon = sequence[index:index + 3]
        aminoAcids.append(CODON_TO_AMINO.get(codon, "?"))
    return "".join(aminoAcids)

# Dictionary with the complementary base of every nucleotide, used by the reverse complement palindromes.
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# reverseComplement: builds the reverse complement of a DNA sequence (the sequence read backwards changing every base for its complementary one)
# sequence: DNA sequence formed by the letters A, C, G and T.
# Returns: the reverse complement sequence.
# Complexity: O(n), where n is the length of the sequence.
def reverseComplement(sequence):
    bases = []
    for index in range(len(sequence) - 1, -1, -1):
        bases.append(COMPLEMENT.get(sequence[index], "N"))
    return "".join(bases)

# longestDnaPalindrome: finds the longest DNA palindrome (reverse complement palindrome) of a sequence, which is a segment equal to its own reverse
# complement, for example GGTACC. These palindromes always have an even length, so the search expands from every position between two bases 
# and compares complementary bases instead of equal ones.
# sequence: DNA sequence to analyze.
# Returns: tuple (start, length) of the longest DNA palindrome; (0, 0) when the sequence does not contain any.
# Complexity: O(n*k) where k is the size of the longest palindrome found. O(n*n) in the worst case.
def longestDnaPalindrome(sequence):
    bestStart = 0
    bestLength = 0
    for center in range(len(sequence) - 1):
        radius = 0
        left = center
        right = center + 1
        while left >= 0 and right < len(sequence) and COMPLEMENT.get(sequence[left], "N") == sequence[right]:
            radius = radius + 1
            left = left - 1
            right = right + 1
        if (2 * radius) > bestLength:
            bestLength = 2 * radius
            bestStart = center - radius + 1
    return bestStart, bestLength