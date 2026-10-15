# Test cases of the algorithms used in the analysis of the SARS-CoV-2 genome.
# Every test compares the result of a function against a value calculated by hand, and at the end the number of passed tests is shown.

import os
import dnaAlgorithms

# Small genome used by the tests, it is read from a file to also test readFasta.
TEST_FILE = "test-genome.txt"

# Counters of passed tests.
passedTests = 0
totalTests = 0

# checkResult: compares an obtained value against the expected one and prints whether the test passed or failed.
# testName: description of the test.
# obtained: value returned by the function that is being tested.
# expected: value calculated by hand.
# Returns: nothing, the counters are updated
# Complexity: O(1).
def checkResult(testName, obtained, expected):
    global passedTests
    global totalTests
    totalTests = totalTests + 1
    if obtained == expected:
        passedTests = passedTests + 1
        print("  PASS  " + testName)
    else:
        print("  FAIL  " + testName + ": obtained " + str(obtained) +
              ", expected " + str(expected))

# testReadFasta: verifies that the header is separated from the sequence and that the sequence lines are joined into a single string.
# Returns: nothing
# Complexity: O(n), with n the size of the test file.
def testReadFasta():
    print("readFasta")
    folder = os.path.dirname(os.path.abspath(__file__))
    header, sequence = dnaAlgorithms.readFasta(os.path.join(folder, "data", TEST_FILE))
    checkResult("the header is the first line", header[0], ">")
    checkResult("the sequence has no line breaks", sequence,
                "ATGGAGAGCCTTACGTAGCATGCA")

# testKmpSearch: verifies the KMP search with patterns that appear once, several times, with overlapping and with patterns that do not appear.
# Returns: nothing
# Complexity: O(n+m) per test.
def testKmpSearch():
    print("kmpSearch")
    genome = "ATGGAGAGCCTTACGTAGCATGCA"
    checkResult("pattern that appears once", dnaAlgorithms.kmpSearch(genome, "AGCCTT"), [6])
    checkResult("pattern that appears twice", dnaAlgorithms.kmpSearch(genome, "GCA"), [17, 21])
    checkResult("pattern that does not appear", dnaAlgorithms.kmpSearch(genome, "GGGG"), [])
    checkResult("overlapped occurrences", dnaAlgorithms.kmpSearch("AAAA", "AA"), [0, 1, 2])
    checkResult("pattern longer than the text", dnaAlgorithms.kmpSearch("ACG", "ACGT"), [])
    checkResult("prefix table of AABAACAABAA",
                dnaAlgorithms.buildPrefixTable("AABAACAABAA"),
                [0, 1, 0, 1, 2, 0, 1, 2, 3, 4, 5])

# testLongestPalindrome: verifies the palindrome search with odd and even
# palindromes, sequences without repetitions and sequences with one letter.
# Returns: nothing
# Complexity: O(n) per test.
def testLongestPalindrome():
    print("longestPalindrome")
    checkResult("odd palindrome CGTGC", dnaAlgorithms.longestPalindrome("GACGTGCTT"), (2, 5))
    checkResult("even palindrome CGGC", dnaAlgorithms.longestPalindrome("GACGGCTT"), (2, 4))
    checkResult("whole sequence is a palindrome",
                dnaAlgorithms.longestPalindrome("AAAA"), (0, 4))
    checkResult("sequence without palindromes longer than one letter",
                dnaAlgorithms.longestPalindrome("ACGT"), (0, 1))
    checkResult("empty sequence", dnaAlgorithms.longestPalindrome(""), (0, 0))

# testLongestCommonSubstring: verifies the dynamic programming solution with sequences that share a substring and with sequences that share nothing.
# Returns: nothingg
# Complexity: O(nm) per test.
def testLongestCommonSubstring():
    print("longestCommonSubstring")
    checkResult("shared substring ACGTA",
                dnaAlgorithms.longestCommonSubstring("ACGTACGT", "TTACGTAA"), (0, 2, 5))
    checkResult("identical sequences",
                dnaAlgorithms.longestCommonSubstring("ACGT", "ACGT"), (0, 0, 4))
    checkResult("sequences without common letters",
                dnaAlgorithms.longestCommonSubstring("AAAA", "CCCC"), (0, 0, 0))

# testTranslateSequence: verifies the translation of DNA into amino acids in the three reading frames and the use of the codon table.
# Returns: nothing
# Complexity: O(n) per test.
def testTranslateSequence():
    print("translateSequence")
    genome = "ATGGAGAGCCTTACGTAGCATGCA"
    checkResult("reading frame 0", dnaAlgorithms.translateSequence(genome, 0), "MESLT*HA")
    checkResult("reading frame 1", dnaAlgorithms.translateSequence(genome, 1), "WRALRSM")
    checkResult("codon of the stop", dnaAlgorithms.CODON_TO_AMINO["TAG"], "*")
    checkResult("codons of the serine", dnaAlgorithms.AMINO_TO_CODONS["S"],
                ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"])

# testDnaPalindrome: verifies the reverse complement and the search of DNA
# palindromes, which are the segments equal to their own reverse complement.
# Returns: nothing
# Complexity: O(n*k) per test.
def testDnaPalindrome():
    print("longestDnaPalindrome")
    checkResult("reverse complement of ATGC", dnaAlgorithms.reverseComplement("ATGC"), "GCAT")
    checkResult("DNA palindrome AAGGTACCTT",
                dnaAlgorithms.longestDnaPalindrome("AAGGTACCTT"), (0, 10))
    checkResult("DNA palindrome ACGT inside the sequence",
                dnaAlgorithms.longestDnaPalindrome("ACGTTT"), (0, 4))
    checkResult("sequence without complementary bases",
                dnaAlgorithms.longestDnaPalindrome("AAAA"), (0, 0))

# main: executes every group of tests and prints the final summary.
# Returns: nothing
# Complexity: dominated by the longest common substring tests, O(n*m).
def main():
    print("TEST CASES")
    testReadFasta()
    testKmpSearch()
    testLongestPalindrome()
    testLongestCommonSubstring()
    testTranslateSequence()
    testDnaPalindrome()
    print("")
    print("Passed tests: " + str(passedTests) + " of " + str(totalTests))

main()
