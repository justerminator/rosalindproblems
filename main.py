from math import comb
import sys, re # IO and motif finding
import urllib3 # get protein fasta from uniprot
from Bio import SeqIO, Seq # fasta parsing and manipulation   


f = open("rosalind_mrna.txt", 'r')

# 7/23/23

def countSymbols(s):

    aCount = 0
    cCount = 0
    gCount = 0
    tCount = 0
    strResult = ""

    for letter in s:
        if letter == 'A':
            aCount += 1
        elif letter == 'C':
            cCount += 1
        elif letter == 'G':
            gCount += 1
        elif letter == 'T':
            tCount += 1

    strResult = str(aCount) + ' ' + str(cCount) + ' ' + str(gCount) + ' ' + str(tCount)

    return strResult

# 7/24/23

# Transcribing DNA into RNA

def dnaTorna(s):

    tReplacement = 'U'
    result = ""

    print(s)
    for letter in s:
        if letter == 'T':
            result += tReplacement
        else:
            result += letter

    return result

# Complementing a Strand of DNA

def revComplementDNA(s):

    aReplacement = 'T'
    tReplacement = 'A'
    gReplacement = 'C'
    cReplacement = 'G'
    result = ""

    for letter in s:
        if letter == 'A':
            result += aReplacement
        elif letter == 'T':
            result += tReplacement
        elif letter == 'G':
            result += gReplacement
        elif letter == 'C':
            result += cReplacement

    return result[::-1]


#Rabbits and Recurrence Relations

def wabbits(month, pairs):

    if month == 1:
        return 1
    if month == 2:
        return 1
    if month < 1:
        return 0


    return (pairs * wabbits(month - 2, pairs)) + wabbits(month - 1, pairs)


# 7/25/23

# Computing GC Content

def loopThrough():

    curr = f.readline()
    numerator = 0
    denominator = 0
    result = 0
    nameStorage = ''
    tempName = ''

    while curr != '' and curr != None:
        if curr[0] == '>':
            tempName = curr
            curr = f.readline()
            while curr != '' and curr[0] != '>':
                numerator += countGC(curr)
                denominator += len(curr.strip())
                #print(numerator)
                #print(denominator)
                curr = f.readline()

            if result < numerator / denominator:
                nameStorage = tempName[1:-1]
                result = numerator / denominator

            numerator = 0
            denominator = 0

        curr = f.readline()

    print(nameStorage)
    print(result * 100)
    return nameStorage, result * 100


def countGC(s) -> int:
    result = 0
    for letter in s:
        if letter == 'G':
            result += 1
        elif letter == 'C':
            result += 1

    return result

# Counting Point Mutations

def hamDistance(s, t):

    result = 0

    for idx in range(len(s)):
        if s[idx] != t[idx]:
            result += 1
    return result

# Mendel's First Law

def mendel(dom, het, rec):

    totPop = dom + het + rec

    totCombos = comb(totPop, 2)
    print(totCombos)

    numDomBabies = comb(dom, 2) + dom*rec + dom*het + .5*het*rec + .75*comb(het, 2)

    prob = numDomBabies / totCombos

    return prob

# Translating RNA into Protein

amino_acids={'AUG':'M','UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S',
                 'UCC':'S','UCA':'S','UCG':'S','CUU': 'L','AUU': 'I','GUU': 'V',
                 'CUC': 'L','AUC': 'I','GUC': 'V','CUA': 'L','AUA': 'I','GUA': 'V',
                 'CUG': 'L','GUG': 'V','CCU': 'P', 'ACU': 'T','GCU': 'A','CCC': 'P',
                 'ACC': 'T','GCC': 'A','CCA': 'P','ACA': 'T','GCA': 'A','CCG': 'P',
                 'ACG': 'T','GCG': 'A','UAU': 'Y','CAU': 'H','AAU': 'N','GAU': 'D',
                 'UAC': 'Y','CAC': 'H','AAC': 'N','GAC': 'D','UAA': '*','CAA': 'Q',
                 'AAA': 'K','GAA': 'E','UAG': '*','CAG': 'Q','AAG': 'K','GAG': 'E',
                 'UGU': 'C','CGU': 'R','AGU': 'S','GGU': 'G','UGC': 'C','CGC': 'R',
                 'AGC': 'S','GGC': 'G','UGA': '*','CGA': 'R','AGA': 'R','GGA': 'G',
                 'UGG': 'W','CGG': 'R','AGG': 'R','GGG': 'G'}

def rnaToProtein(s):
    lenString = len(s) - 2

    result = ""

    i = 0
    while i < len(s):
        #print(s[i:i + 3])
        if s[i:i + 3] in amino_acids:

            result += amino_acids[s[i:i + 3]]
            i += 2
        i += 1

    return result

# 7/26/23

# Finding a Motif in DNA

def findMotif(s, t) -> list:
    result = []
    tRange = len(t)
    toPrint = ""
    for i in range(len(s)):
        if s[i] == t[0]:
            if s[i:i + tRange] == t:
                result.append(i + 1)
                #toPrint += str(i + 1) + " "
    #print(toPrint)
    return result


# 7/27/23

# Consensus and Profile

def conAndpro():
    curr = f.readline()
    curr = f.readline()

    temp = curr
    length = 0
    stopper = 0
    entireSequence = ""

    # while temp != '' and temp != None:
    #     if temp[0] == '>':
    #         if stopper == 0:
    #             stopper += 1
    #         else:
    #             break
    #         temp = f.readline()
    #     for idx in range(len(temp)):
    #         print(temp)
    #         length += len(temp.strip())
    #         temp = f.readline()
    #         if temp[0] == '>':
    #             break
    #     break
    #
    #
    # print(length)

    result = ""

    stopper = 0
    tracker = 0
    stop2 = False
    while curr != '' and curr != None:
        if curr[0] == '>':
            if stopper == 0 and stop2:
                stopper += 1
            elif stop2 == False:
                tracker = len(entireSequence)
                stop2 = True
            curr = f.readline()
        for letter in curr.strip():
            entireSequence += str(letter)
        print(curr)
        curr = f.readline()

    print("entire: ", entireSequence)
    print("tracker ", tracker)
    a = [0] * tracker
    t = [0] * tracker
    c = [0] * tracker
    g = [0] * tracker

    lineCounter = 0
    for idx in range(len(entireSequence)):
        if idx % tracker == 0:
            lineCounter = 0
        if entireSequence[idx] == 'A':
            a[lineCounter] += 1
        if entireSequence[idx] == 'T':

            t[lineCounter] += 1
        if entireSequence[idx] == 'C':
            c[lineCounter] += 1
        if entireSequence[idx] == 'G':
            print("bruh ", lineCounter)
            g[lineCounter] += 1


        lineCounter += 1





    for x in range(len(a)):
        if a[x] >= t[x] and a[x] >= c[x] and a[x] >= g[x]:
            result += "A"
        elif t[x] >= a[x] and t[x] >= c[x] and t[x] >= g[x]:
            result += "T"
        elif g[x] >= a[x] and g[x] >= c[x] and g[x] >= t[x]:
            result += "G"
        elif c[x] >= a[x] and c[x] >= g[x] and c[x] >= t[x]:
            result += "C"

    print(result)

    print("A:", *a)
    print("C:", *c)
    print("G:", *g)
    print("T:", *t)

    print(len(a), len(t), len(c), len(g))

# 7/28/23

# Mortal Fibonacci Rabbits

def mortalWabbits(month, lifeSpan):

    sequence = [1, 1]
    for i in range(month - 2):
        newNum = 0
        if i + 2 < lifeSpan:
            newNum = sequence[i] + sequence[i + 1]

        else:
            for j in range(lifeSpan - 1):
                newNum += sequence[i - j]
                #print(sequence[i - j])
            #print('bruh')


        sequence.append(newNum)

    print(sequence)



# 7/29/23

# Overlap Graphs

def dataAcquirer():
    nameNodes = []
    listS = []
    toAdd = ""

    curr = f.readline()
    while curr != '' and curr != None:
        if curr[0] == '>':
            nameNodes.append(curr[1:].strip())
            if toAdd != "":
                listS.append(toAdd)
                toAdd = ""

        else:
            toAdd += curr.strip()
        curr = f.readline()
    listS.append(toAdd)

    #subFinder(nameNodes, listS, 3)

    return nameNodes, listS

def subFinder(nameNodes, listS, k):
    #print(listS)
    #print(nameNodes)

    lenListS = len(listS)
    result = []

    for i in range(lenListS):
        toCompare = listS[i][len(listS[i]) - k:]
        #print(listS[i][len(listS[i]) - k:])
        for j in range(i + 1, lenListS):
            if toCompare == listS[j][:3]:
                if (nameNodes[i], nameNodes[j]) not in result:
                    result.append((nameNodes[i], nameNodes[j]))
        for x in range(0, i):
            if toCompare == listS[x][:3]:
                if (nameNodes[i], nameNodes[x]) not in result:
                    result.append((nameNodes[i], nameNodes[x]))

    # print result
    for i in range(len(result)):
        print(result[i][0], result[i][1])

# Calculating Expected Offspring

def expectedOffspring() -> float:
    curr = f.readline()
    dataset = []

    toAdd = ''
    for letter in curr.strip():
        if letter != ' ':
            toAdd += letter
        elif letter == ' ':
            dataset.append(int(toAdd))
            toAdd = ''
    dataset.append(int(toAdd))

    #print(dataset)
    result = 2 * dataset[0] + 2 * dataset[1] + 2 * dataset[2] + 2 * .75*dataset[3] + 2 * 0.5*dataset[4]
    #print(result)

    return result

# Finding a Shared Motif

def sharedMotif():
    nameSequences, listS = dataAcquirer()

    listS.sort()
    #print(listS)

    toCompare = listS[0]

    motif = ''

    for i in range(len(toCompare)):
        n = 0
        present = True
        while present:
            for idx in range(1, len(listS)):
                if toCompare[i:i + n] not in listS[idx] or n > len(toCompare):
                    present = False
                    break

            if present:
                motif = max(toCompare[i:i+n], motif, key=len)
            n += 1
    return motif



# 7/31/23

# Independent Alleles


def indAlleles(k, n):

    totPop = 2 ** k

    result = 0
    for i in range(n, (totPop + 1)):
        result += comb(totPop, i) * 0.25**i * 0.75**(totPop-i)
    
    return result


# 8/1/23

# Inferring mRNA from Protein

def RNAtoAA():
    frequence = {}
    for k, v in amino_acids.items():
        if v not in frequence:
            frequence[v] = 0
        frequence[v] += 1
    
    return frequence

def possiblemRNA():
    freq = RNAtoAA()
    possibilityAmnt = 1
    print(freq)
    for letter in f.readline().strip():
        
        possibilityAmnt = (possibilityAmnt * freq[letter]) % 1000000
    
    possibilityAmnt = (possibilityAmnt * freq['*']) % 1000000

    return possibilityAmnt

print(possiblemRNA())

f.close()