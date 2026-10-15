# SARS-CoV-2 genome analysis.
# Authors: Felipe Merida - A00838419, Isaac Pereyra - A00841493, Diego Hanger - A00843443
# Date: 2026-10-14

import os
import dnaAlgorithms

WUHAN_FILE = "SARS-COV-2-MN908947.3.txt"
TEXAS_FILE = "SARS-COV-2-MT106054.1.txt"
PROTEINS_FILE = "seq-proteins.txt"

GENE_FILES = [
    ("M", "gen-M.txt"),
    ("S", "gen-S.txt"),
    ("ORF1ab", "gen-ORF1AB.txt")
]

PROTEIN_GENES = {
    "QHD43416": "S",
    "QHD43417": "ORF3a",
    "QHD43418": "E",
    "QHD43419": "M",
    "QHD43420": "ORF6",
    "QHD43421": "ORF7a",
    "QHD43422": "ORF8",
    "QHD43423": "N",
    "QHI42199": "ORF10"
}


# Here we get the path for a file in the data folder.
# fileName: name of the file.
# Returns: the path to that file.
# Time: O(1) for these fixed paths, Space: O(1).
def buildPath(fileName):
    folder = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(folder, "data", fileName)


# First we look for each gene using KMP.
# genome: Wuhan DNA sequence.
# genes: list of gene names and sequences
# Returns: list with each gene name and its positions.
# Time: O(g*n + total gene lengths), Space: O(m + r) per search.
def findGenes(genome, genes):
    located = []

    print("Genes")

    for name, sequence in genes:
        positions = dnaAlgorithms.kmpSearch(genome, sequence)

        if len(positions) == 0:
            print(name, "not found")
        else:
            start = positions[0]
            end = start + len(sequence) - 1

            located.append((name, start, end))

            print(name, "start:", start, "end:", end)

    return located


# Finds the longest palindrome in each gene.
# genes: gene names and their DNA sequences.
# located: positions of the genes in the genome.
# Returns: nothing, prints the palindromes and positions.
# Time: O(sum of gene lengths squared), Space: O(m) for printing.
def findPalindromes(genes, located):
    print("Palindromes")

    for name, sequence in genes:
        geneStart = -1

        for geneName, start, _ in located:
            if geneName == name:
                geneStart = start

        start, length = dnaAlgorithms.longestPalindrome(sequence)
        print(name, sequence[start:start + length], "start:", start, "length:", length)

        if geneStart >= 0:
            print("Genome position:", geneStart + start)


# So, this tries a protein that changes its reading frame once.
# genome: DNA sequence.
# protein: amino acid sequence.
# translations: DNA translated in the three frames.
# Returns: (start, end, shiftAfter) or None.
# Time: O(n + c*p^2) worst case, Space: O(c).
def findShiftedProtein(genome, protein, translations):
    prefixSize = min(8, len(protein))
    prefix = protein[:prefixSize]

    for frame in range(3):
        positions = dnaAlgorithms.kmpSearch(translations[frame], prefix)

        for position in positions:
            start = frame + 3 * position

            # How many amino acids match without changing the frame.
            normalMatches = 0
            matching = True

            while matching and normalMatches < len(protein):
                codonStart = start + 3 * normalMatches
                codon = genome[codonStart:codonStart + 3]

                if dnaAlgorithms.CODON_TO_AMINO.get(codon) == protein[normalMatches]:
                    normalMatches += 1
                else:
                    matching = False

            # We go back one base after each amino acid and check the rest.
            for shiftAfter in range(1, normalMatches + 1):
                total = shiftAfter
                matching = True

                while matching and total < len(protein):
                    codonStart = start + 3 * total - 1
                    codon = genome[codonStart:codonStart + 3]

                    if dnaAlgorithms.CODON_TO_AMINO.get(codon) == protein[total]:
                        total += 1
                    else:
                        matching = False

                if total == len(protein):
                    end = start + 3 * len(protein) - 2
                    return start, end, shiftAfter

    return None


# Here we check where each protein is encoded.
# genome: Wuhan DNA sequence.
# proteins: protein names and amino acid sequences.
# Returns: list of protein names, positions and frame changes.
# Time: O(k*n + sum(c*p^2)), Space: O(n + k).
def findProteins(genome, proteins):
    print("Proteins")

    translations = []

    for frame in range(3):
        translations.append(dnaAlgorithms.translateSequence(genome, frame))

    regions = []

    for name, protein in proteins:
        result = None

        for frame in range(3):
            positions = dnaAlgorithms.kmpSearch(translations[frame], protein)

            if len(positions) > 0 and result is None:
                start = frame + 3 * positions[0]
                end = start + 3 * len(protein) - 1
                result = (start, end, -1)

        if result is None:
            result = findShiftedProtein(genome, protein, translations)

        if result is None:
            print(name, "not found")
        else:
            start, end, shiftAfter = result

            regions.append((name, start, end, shiftAfter))

            print(name, "start:", start, "end:", end)

            if shiftAfter >= 0:
                print(name, "frame shift after amino acid:", shiftAfter)

    return regions


# Checks which protein contains a given position.
# position: nucleotide index.
# regions: protein regions found previously.
# Returns: protein region or None
# Time: O(k), Space: O(1).
def findProteinOfPosition(position, regions):
    for region in regions:
        if region[1] <= position <= region[2]:
            return region

    return None


# Gets the gene name, based on the protein ID.
# proteinName: protein identifier.
# Returns: gene name.
# Time: O(1), Space: O(1).
def findGeneName(proteinName):
    if proteinName[:8] == "QHD43415":
        return "ORF1ab"

    return PROTEIN_GENES.get(proteinName, proteinName)


# We need the first base of the codon, not the changed base.
# position: index of the changed base.
# region: protein start, end and reading-frame information.
# Returns: starting position of its codon.
# Time: O(1), Space: O(1).
def getCodonStart(position, region):
    start = region[1]
    shiftAfter = region[3]

    if shiftAfter >= 0 and position >= start + 3 * shiftAfter - 1:
        firstShifted = start + 3 * shiftAfter - 1
        return firstShifted + 3 * ((position - firstShifted) // 3)

    return start + 3 * ((position - start) // 3)


# Basically, we compare both genomes and the amino acids.
# wuhan: original Wuhan DNA sequence.
# texas: Texas DNA sequence.
# regions: identified protein regions.
# Returns: list with the differences and their effects.
# Time: O(n + d*k), Space: O(d).
def compareGenomes(wuhan, texas, regions):
    print("Genome comparison")

    size = min(len(wuhan), len(texas))
    differences = []
    equalStart = 0
    results = []

    for i in range(size):
        if wuhan[i] != texas[i]:
            differences.append(i)

    print("Wuhan length:", len(wuhan))
    print("Texas length:", len(texas))
    print("Equal bases:", size - len(differences))
    print("Different bases:", len(differences))

    for position in differences:
        if equalStart < position:
            print("Equal region:", equalStart, position - 1)

        equalStart = position + 1

    if equalStart < size:
        print("Equal region:", equalStart, size - 1)

    for position in differences:
        region = findProteinOfPosition(position, regions)

        print("Difference:", position, wuhan[position], texas[position])

        if region is None:
            print("Outside identified coding regions")
            results.append((position, "outside", "not classified"))

        else:
            gene = findGeneName(region[0])
            codonStart = getCodonStart(position, region)

            wuhanCodon = wuhan[codonStart:codonStart + 3]
            texasCodon = texas[codonStart:codonStart + 3]

            commonStart, _, commonLength = dnaAlgorithms.longestCommonSubstring(
                wuhanCodon, texasCodon)
            commonPart = wuhanCodon[commonStart:commonStart + commonLength]

            wuhanAmino = dnaAlgorithms.CODON_TO_AMINO.get(wuhanCodon, "?")
            texasAmino = dnaAlgorithms.CODON_TO_AMINO.get(texasCodon, "?")

            if wuhanAmino == "?" or texasAmino == "?":
                effect = "unknown"
            elif wuhanAmino == texasAmino:
                effect = "silent"
            elif wuhanAmino == "*" or texasAmino == "*":
                effect = "stop codon changed"
            else:
                effect = "missense"

            print("Gene:", gene, "codons:", wuhanCodon, texasCodon,
                  "amino acids:", wuhanAmino, texasAmino, "effect:", effect)
            print("Longest common substring in codons:", commonPart)

            results.append((position, gene, effect))

    if len(wuhan) != len(texas):
        print("Length difference:", abs(len(wuhan) - len(texas)))
        print("Remaining bases not compared")

    return results


# Reads the files and runs the four parts of the project.
# Parameters: none.
# Returns: nothing, prints results.
# Time: dominated by palindrome searches and frame-shift search.
# Space: O(n) for the genomes, genes and proteins in memory.
def main():
    _, wuhan = dnaAlgorithms.readFasta(buildPath(WUHAN_FILE))
    _, texas = dnaAlgorithms.readFasta(buildPath(TEXAS_FILE))

    genes = []

    for name, fileName in GENE_FILES:
        _, sequence = dnaAlgorithms.readFasta(buildPath(fileName))
        genes.append((name, sequence))

    proteins = dnaAlgorithms.readFastaMultiple(buildPath(PROTEINS_FILE))

    print("Positions start at 0")

    located = findGenes(wuhan, genes)
    findPalindromes(genes, located)
    regions = findProteins(wuhan, proteins)
    compareGenomes(wuhan, texas, regions)


if __name__ == "__main__":
    main()