# Main program of the integrated activity: analysis of the SARS-CoV-2 genome.
# It locates the M, S and ORF1ab genes inside the Wuhan genome, finds the longest palindrome of each gene, locates the regions that encode 
# the supplied proteins and compares the Wuhan and Texas genomes.
import os
import dnaAlgorithms

# Names of the supplied input files; they must be in the same folder as this file.
WUHAN_FILE = "SARS-COV-2-MN908947.3.txt"
TEXAS_FILE = "SARS-COV-2-MT106054.1.txt"
PROTEINS_FILE = "seq-proteins.txt"
GENE_FILES = [("M", "gen-M.txt"), ("S", "gen-S.txt"), ("ORF1AB", "gen-ORF1AB.txt")]

# buildPath: builds the absolute path of an input file, so the program can be executed from any folder.
# fileName: name of the file that is going to be read, stored in the folder
# "data" next to this program.
# Returns: absolute path of the file.
# Complexity: O(1).
def buildPath(fileName):
    folder = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(folder, "data", fileName)

# findGenes: locates every supplied gene inside the genome with the KMP algorithm and prints the positions where each gene starts and ends.
# genome: complete genome sequence.
# genes: list of tuples (geneName, geneSequence).
# Returns: list of tuples (geneName, start, end) with 0 based positions.
# Complexity: O(k*(n+m)) with k genes, genome of size n and gene of size m.
def findGenes(genome, genes):
    located = []
    print("1) GENE IDENTIFICATION (Wuhan 2019 genome)")
    for geneName, geneSequence in genes:
        positions = dnaAlgorithms.kmpSearch(genome, geneSequence)
        if len(positions) == 0:
            print("  Gene " + geneName + ": not found in the genome")
        else:
            start = positions[0]
            end = start + len(geneSequence) - 1
            located.append((geneName, start, end))
            print("  Gene " + geneName + ": start " + str(start) + ", end " + str(end) +
                  ", length " + str(len(geneSequence)) + ", occurrences " + str(len(positions)))
            print("    first 12 bases: " + geneSequence[0:12])
    print("")
    return located

# findPalindromes: finds and prints the longest palindrome of every gene.
# genes: list of tuples (geneName, geneSequence).
# located: list of tuples (geneName, start, end) returned by findGenes, used to report the position of the palindrome inside the complete genome.
# Returns: nothing, the results are printed.
# Complexity: O(k*m) with k genes and m the length of each gene.
def findPalindromes(genes, located):
    print("2) LONGEST PALINDROME OF EACH GENE")
    for geneName, geneSequence in genes:
        start, length = dnaAlgorithms.longestPalindrome(geneSequence)
        genomeStart = -1
        for name, geneStart, geneEnd in located:
            if name == geneName:
                genomeStart = geneStart + start
        print("  Gene " + geneName + ": length " + str(length) +
              ", position in the gene " + str(start) + " - " + str(start + length - 1) +
              ", position in the genome " + str(genomeStart) + " - " +
              str(genomeStart + length - 1))
        print("    " + geneSequence[start:start + length])
        # A DNA palindrome is a segment equal to its own reverse complement,
        # which is the kind of palindrome that has a biological meaning.
        dnaStart, dnaLength = dnaAlgorithms.longestDnaPalindrome(geneSequence)
        dnaGenomeStart = -1
        if genomeStart >= 0:
            dnaGenomeStart = (genomeStart - start) + dnaStart
        print("    longest DNA palindrome (reverse complement): length " + str(dnaLength) +
              ", position in the gene " + str(dnaStart) + " - " + str(dnaStart + dnaLength - 1) +
              ", position in the genome " + str(dnaGenomeStart) + " - " +
              str(dnaGenomeStart + dnaLength - 1))
        print("    " + geneSequence[dnaStart:dnaStart + dnaLength])
    print("")

# findProteins: determines which region of the genome encodes every supplied protein. The genome is translated in the three reading frames
# and the protein is searched inside the translation with the KMP algorithm.
# genome: complete genome sequence.
# proteins: list of tuples (proteinName, aminoAcidSequence).
# Returns: nothing
# Complexity: O(n+k*(n/3+p)) with n the genome length, k the number of proteins and p the length of each protein.
def findProteins(genome, proteins):
    print("3) REGIONS THAT ENCODE EACH PROTEIN")
    translations = []
    for frame in range(3):
        translations.append(dnaAlgorithms.translateSequence(genome, frame))
    for proteinName, proteinSequence in proteins:
        frame = 0
        start = -1
        while frame < 3 and start < 0:
            positions = dnaAlgorithms.kmpSearch(translations[frame], proteinSequence)
            if len(positions) > 0:
                start = frame + (3 * positions[0])
            frame = frame + 1
        if start < 0:
            print("  " + proteinName + ": not found, no reading frame of the genome " +
                  "encodes this protein")
        else:
            end = start + (3 * len(proteinSequence)) - 1
            print("  " + proteinName + ": genome " + str(start) + " - " + str(end) +
                  ", reading frame " + str(start % 3) + ", " + str(len(proteinSequence)) +
                  " amino acids")
            print("    first 12 bases: " + genome[start:start + 12] +
                  " -> " + proteinSequence[0:4])
    print("")

# findGeneOfPosition: tells which of the located genes contains a position.
# position: index inside the genome.
# located: list of tuples (geneName, start, end).
# Returns: name of the gene, or "outside the supplied genes" when the position does not belong to any of them.
# Complexity: O(k), with k the number of located genes.
def findGeneOfPosition(position, located):
    geneName = "outside the supplied genes"
    for name, start, end in located:
        if position >= start and position <= end:
            geneName = "gene " + name
    return geneName

# compareGenomes: compares the two genomes position by position, reports the regions where they are identical, the positions where 
#they differ and whether each difference changes the encoded amino acid.
# wuhan: genome of Wuhan 2019.
# texas: genome of Texas 2020.
# located: list of tuples (geneName, start, end) of the genes in the genome.
# Returns: nothing, the results are printed.
# Complexity: O(n), with n the length of the shortest genome.
def compareGenomes(wuhan, texas, located):
    print("4) COMPARISON BETWEEN THE WUHAN 2019 AND TEXAS 2020 GENOMES")
    shortest = min(len(wuhan), len(texas))
    differences = []
    for index in range(shortest):
        if wuhan[index] != texas[index]:
            differences.append(index)
    print("  Length of the Wuhan genome: " + str(len(wuhan)))
    print("  Length of the Texas genome: " + str(len(texas)))
    print("  Compared positions: " + str(shortest) + ", equal: " +
          str(shortest - len(differences)) + ", different: " + str(len(differences)))

    print("  Regions where both genomes are identical:")
    regionStart = 0
    for index in differences:
        if index > regionStart:
            print("    " + str(regionStart) + " - " + str(index - 1) +
                  " (" + str(index - regionStart) + " bases)")
        regionStart = index + 1
    if regionStart < shortest:
        print("    " + str(regionStart) + " - " + str(shortest - 1) +
              " (" + str(shortest - regionStart) + " bases)")

    # For every difference the affected codon is the one that starts at th position that changed, which is the criterion used in the expected
    # results of the activity.
    print("  Positions where the genomes differ:")
    for index in differences:
        wuhanCodon = wuhan[index:index + 3]
        texasCodon = texas[index:index + 3]
        wuhanAmino = dnaAlgorithms.CODON_TO_AMINO.get(wuhanCodon, "?")
        texasAmino = dnaAlgorithms.CODON_TO_AMINO.get(texasCodon, "?")
        effect = "the amino acid changes (missense mutation)"
        if wuhanAmino == texasAmino:
            effect = "the amino acid does not change (silent mutation)"
        elif wuhanAmino == "*" or texasAmino == "*":
            effect = "a stop codon appears or disappears (nonsense mutation)"
        print("    index " + str(index) + " (" + findGeneOfPosition(index, located) + "): " +
              "Wuhan " + wuhan[index] + " vs Texas " + texas[index])
        print("      codon " + wuhanCodon + " (" + wuhanAmino + ") vs " +
              texasCodon + " (" + texasAmino + "): " + effect)
    if len(wuhan) != len(texas):
        print("    the Wuhan genome has " + str(len(wuhan) - shortest) +
              " extra bases at the end, positions " + str(shortest) + " - " +
              str(len(wuhan) - 1))

    print("")

# main: reads the input files and executes the four functionalities.
# Returns: nothing
# Complexity: dominated by the palindrome search
def main():
    wuhanHeader, wuhanGenome = dnaAlgorithms.readFasta(buildPath(WUHAN_FILE))
    texasHeader, texasGenome = dnaAlgorithms.readFasta(buildPath(TEXAS_FILE))
    genes = []
    for geneName, geneFile in GENE_FILES:
        geneHeader, geneSequence = dnaAlgorithms.readFasta(buildPath(geneFile))
        genes.append((geneName, geneSequence))
    proteins = dnaAlgorithms.readFastaMultiple(buildPath(PROTEINS_FILE))

    print("SARS-CoV-2 GENOME ANALYSIS")
    print("Wuhan 2019 file: " + WUHAN_FILE + " (" + str(len(wuhanGenome)) + " bases)")
    print("Texas 2020 file: " + TEXAS_FILE + " (" + str(len(texasGenome)) + " bases)")
    print("All the positions are 0 based; the biological numbering adds 1 to them.")
    print("")

    located = findGenes(wuhanGenome, genes)
    findPalindromes(genes, located)
    findProteins(wuhanGenome, proteins)
    compareGenomes(wuhanGenome, texasGenome, located)

main()
