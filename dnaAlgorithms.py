# Algorithms for the SARS-CoV-2 genome analysis project
# Authors: Felipe Merida - A00838419, Isaac Pereyra - A00841493, Diego Hanger - A00843443
# Last modified: 2026-10-14


# Here are the codons that can produce each amino acid.
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


# So this builds the table to get an amino acid from a codon.
# Parameters: none.
# Returns: a dictionary with codons as keys and amino acids as values.
# Time: O(1), Space: O(1), because there are 64 codons.
def buildCodonToAmino():
    table = {}
    for amino in AMINO_TO_CODONS:
        for codon in AMINO_TO_CODONS[amino]:
            table[codon] = amino
    return table


CODON_TO_AMINO = buildCodonToAmino()


# Reads one FASTA record, joining the DNA lines together.
# filePath: path of the file to read.
# Returns: the header and the DNA sequence without line breaks.
# Time: O(n), Space: O(n), where n is the file size.
def readFasta(filePath):
    header = ""
    pieces = []

    with open(filePath, "r") as inputFile:
        header = inputFile.readline().strip()
        for line in inputFile:
            piece = line.strip()
            if piece != "":
                pieces.append(piece)

    sequence = "".join(pieces)
    return header, sequence


# Similar to readFasta, but here the file has several records.
# filePath: path of the file to read.
# Returns: a list of (name, sequence) pairs.
# Time: O(n), Space: O(n), where n is the file size.
def readFastaMultiple(filePath):
    records = []
    currentName = ""
    currentPieces = []

    with open(filePath, "r") as inputFile:
        for line in inputFile:
            cleanLine = line.strip()
            if cleanLine != "":
                if cleanLine[0] == ">":
                    if currentName != "":
                        records.append((currentName, "".join(currentPieces)))
                    currentName = cleanLine[1:]
                    currentPieces = []
                elif currentName != "":
                    currentPieces.append(cleanLine)

    if currentName != "":
        records.append((currentName, "".join(currentPieces)))

    return records


# KMP uses this table so we do not start matching from zero again.
# pattern: sequence used to calculate the matching prefixes.
# Returns: a list with the prefix length for each position.
# Time: O(m), Space: O(m), where m is the pattern length
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
            index = index + 1

    return table


# Searches the text with KMP, including repeated matches
# text: sequence where the search is done.
# pattern: sequence to find.
# Returns: a list of starting positions, using 0-based indexes.
# Time: O(n + m), Space: O(m + r), where r is the number of matches.
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


# Checks palindromes from each center, for odd and even lengths.
# sequence: text or DNA sequence to check.
# Returns: (start, length) of the longest palindrome
# Time: O(n^2) in the worst case, Space: O(1).
def longestPalindrome(sequence):
    bestStart = 0
    bestLength = 0

    for center in range(len(sequence)):
        left = center
        right = center

        while left >= 0 and right < len(sequence) and sequence[left] == sequence[right]:
            left = left - 1
            right = right + 1

        if right - left - 1 > bestLength:
            bestLength = right - left - 1
            bestStart = left + 1

        left = center
        right = center + 1

        while left >= 0 and right < len(sequence) and sequence[left] == sequence[right]:
            left = left - 1
            right = right + 1

        if right - left - 1 > bestLength:
            bestLength = right - left - 1
            bestStart = left + 1

    return bestStart, bestLength


# Here we use DP to find the longest shared part of two strings.
# first: first sequence to compare.
# second: second sequence to compare.
# Returns: (start in first, start in second, length).
# Time: O(n * m), Space: O(m), using two DP rows.
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


# Basically we read three DNA bases at a time to get each codon.
# sequence: DNA sequence to translate.
# frame: reading frame (0, 1 or 2)
# Returns: the amino acid sequence, with * for stop codons.
# Time: O(n), Space: O(n), for the translated sequence.
def translateSequence(sequence, frame):
    aminoAcids = []

    for index in range(frame, len(sequence) - 2, 3):
        codon = sequence[index:index + 3]
        aminoAcids.append(CODON_TO_AMINO.get(codon, "?"))

    return "".join(aminoAcids)