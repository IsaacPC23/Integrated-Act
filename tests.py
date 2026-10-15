# Tests for the SARS-CoV-2 genome analysis project.
# Authors: Felipe Merida - A00838419, Isaac Pereyra - A00841493, Diego Hanger - A00843443
# Date: 2026-10-14

import dnaAlgorithms
import main

passedTests = 0
totalTests = 0


# So, this checks if a result is the same as what we expected.
# testName: name of the test.
# obtained: result from the function.
# expected: result we want to get.
# Returns: nothing, updates the counters.
# Time: O(n) for comparing sequences, Space: O(1) extra if they match.
def checkResult(testName, obtained, expected):
    global passedTests
    global totalTests

    totalTests += 1

    if obtained == expected:
        passedTests += 1
        print("PASS", testName)
    else:
        print("FAIL", testName, "got:", obtained, "expected:", expected)


# Tests if FASTA lines are joined into one sequence.
# Parameters: none.
# Returns: nothing, checks the header and sequence.
# Time: O(n), Space: O(n), where n is the file size.
def testReadFasta():
    header, sequence = dnaAlgorithms.readFasta(main.buildPath("test-genome.txt"))

    checkResult("FASTA header", header[0], ">")
    checkResult("FASTA sequence", sequence, "ATGGAGAGCCTTACGTAGCATGCA")


# Here we test KMP with some different cases.
# Parameters: none.
# Returns: nothing, checks the search results.
# Time: O(n + m) for each search, Space: O(m + r).
def testKmpSearch():
    genome = "ATGGAGAGCCTTACGTAGCATGCA"

    checkResult("one match", dnaAlgorithms.kmpSearch(genome, "AGCCTT"), [6])
    checkResult("two matches", dnaAlgorithms.kmpSearch(genome, "GCA"), [17, 21])
    checkResult("no match", dnaAlgorithms.kmpSearch(genome, "GGGG"), [])
    checkResult("overlapping matches", dnaAlgorithms.kmpSearch("AAAA", "AA"), [0, 1, 2])
    checkResult("pattern too long", dnaAlgorithms.kmpSearch("ACG", "ACGT"), [])
    checkResult("prefix table", dnaAlgorithms.buildPrefixTable("AABAACAABAA"),
                [0, 1, 0, 1, 2, 0, 1, 2, 3, 4, 5])


# Tests palindromes with odd, even and empty sequences.
# Parameters: none.
# Returns: nothing, checks the palindrome positions.
# Time: O(n^2) worst case, Space: O(1) extra.
def testLongestPalindrome():
    checkResult("odd palindrome", dnaAlgorithms.longestPalindrome("GACGTGCTT"), (2, 5))
    checkResult("even palindrome", dnaAlgorithms.longestPalindrome("GACGGCTT"), (2, 4))
    checkResult("whole palindrome", dnaAlgorithms.longestPalindrome("AAAA"), (0, 4))
    checkResult("one letter palindrome", dnaAlgorithms.longestPalindrome("ACGT"), (0, 1))
    checkResult("empty palindrome", dnaAlgorithms.longestPalindrome(""), (0, 0))


# Basically, this checks the shared substring using DP.
# Parameters: none.
# Returns: nothing, checks the common substring.
# Time: O(n*m), Space: O(m).
def testLongestCommonSubstring():
    checkResult("shared substring", dnaAlgorithms.longestCommonSubstring(
        "ACGTACGT", "TTACGTAA"), (0, 2, 5))

    checkResult("same sequences", dnaAlgorithms.longestCommonSubstring(
        "ACGT", "ACGT"), (0, 0, 4))

    checkResult("no shared substring", dnaAlgorithms.longestCommonSubstring(
        "AAAA", "CCCC"), (0, 0, 0))


# This checks DNA translation in two frames and some codons.
# Parameters: none.
# Returns: nothing, checks amino acid values.
# Time: O(n), Space: O(n) per translation.
def testTranslateSequence():
    genome = "ATGGAGAGCCTTACGTAGCATGCA"

    checkResult("frame 0", dnaAlgorithms.translateSequence(genome, 0), "MESLT*HA")
    checkResult("frame 1", dnaAlgorithms.translateSequence(genome, 1), "WRALRSM")
    checkResult("stop codon", dnaAlgorithms.CODON_TO_AMINO["TAG"], "*")
    checkResult("serine codons", dnaAlgorithms.AMINO_TO_CODONS["S"],
                ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"])


# Checks if we find a gene and what happens if it is missing.
# Parameters: none.
# Returns: nothing, checks gene start and end positions.
# Time: O(n + m) per gene, Space: O(m + r).
def testFindGenes():
    found = main.findGenes("GGGATGCCCTTT", [("smallGene", "ATGCCC")])
    checkResult("gene found", found, [("smallGene", 3, 8)])

    missing = main.findGenes("ACGT", [("smallGene", "TTTT")])
    checkResult("gene missing", missing, [])


# Here we check a small protein and the real ORF1ab example.
# Parameters: none.
# Returns: nothing, checks protein regions and a frame shift.
# Time: depends on KMP and the shifted protein search, Space: O(n + k).
def testFindProteins():
    found = main.findProteins("CCCATGGAAAGCTAA", [("smallProtein", "MES")])

    checkResult("small protein", found, [("smallProtein", 3, 11, -1)])

    _, wuhan = dnaAlgorithms.readFasta(main.buildPath(main.WUHAN_FILE))
    proteins = dnaAlgorithms.readFastaMultiple(main.buildPath(main.PROTEINS_FILE))
    regions = main.findProteins(wuhan, proteins)

    checkResult("24 proteins found", len(regions), 24)

    shifted = None

    for region in regions:
        if region[0] == "QHD43415_11":
            shifted = region

    checkResult("protein with frame shift", shifted,
                ("QHD43415_11", 13441, 16235, 9))


# Checks silent and missense mutations using the right codons.
# Parameters: none.
# Returns: nothing, checks mutation positions and effects.
# Time: O(n + d*k), Space: O(d).
def testCompareGenomes():
    wuhan = "ATGGCTTTC"
    texas = "ATGGTTTTT"
    regions = [("smallProtein", 0, 8, -1)]

    result = main.compareGenomes(wuhan, texas, regions)

    checkResult("small genome mutations", result,
                [(4, "smallProtein", "missense"),
                 (8, "smallProtein", "silent")])

    _, wuhan = dnaAlgorithms.readFasta(main.buildPath(main.WUHAN_FILE))
    _, texas = dnaAlgorithms.readFasta(main.buildPath(main.TEXAS_FILE))

    proteins = dnaAlgorithms.readFastaMultiple(main.buildPath(main.PROTEINS_FILE))
    regions = main.findProteins(wuhan, proteins)
    result = main.compareGenomes(wuhan, texas, regions)

    positions = []
    silent = 0
    missense = 0

    for position, gene, effect in result:
        positions.append(position)

        if effect == "silent":
            silent += 1
        elif effect == "missense":
            missense += 1

    checkResult("real genome positions", positions,
                [8781, 18602, 18974, 19174, 27924, 28143, 29094])

    checkResult("silent mutations", silent, 4)
    checkResult("missense mutations", missense, 3)


# Runs all the tests and shows how many passed.
# Parameters: none.
# Returns: nothing, prints the final count.
# Time: depends on the most expensive tests, Space: O(n + k).
def mainTests():
    testReadFasta()
    testKmpSearch()
    testLongestPalindrome()
    testLongestCommonSubstring()
    testTranslateSequence()
    testFindGenes()
    testFindProteins()
    testCompareGenomes()

    print("Passed tests:", passedTests, "of", totalTests)

    if passedTests != totalTests:
        raise AssertionError("Some tests failed")


if __name__ == "__main__":
    mainTests()